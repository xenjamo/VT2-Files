#include "NEOM9N_thread.h"
#include <cstdint>

NEOM9N::NEOM9N(PinName Tx, PinName Rx) :
    m_bufferedSerial(Tx, Rx),
    thread(osPriorityNormal, 4096)
{
    m_bufferedSerial.set_baud(38400);
    m_bufferedSerial.set_format(8, BufferedSerial::None, 1);
    m_bufferedSerial.set_blocking(false);
#if PRINT_FOR_DEBUG
    m_run_timer.start();
#endif
    m_pos_ecef_0.setZero();
    m_pos_ecef.setZero();
    m_pos_enu.setZero();
    m_R_ecefToLocal_0.setIdentity();
}

NEOM9N::~NEOM9N()
{
    ticker.detach();
}

void NEOM9N::start_loop()
{
    thread.start(callback(this, &NEOM9N::update));
    ticker.attach(callback(this, &NEOM9N::sendThreadFlag), std::chrono::microseconds{ static_cast<long int>( 1.0e6f * m_Ts ) });
}

void NEOM9N::reset_local()
{
    // copy so it can not be updated during calculations (maybe mutex this)
    ubxNavPVT_t ubxNavPVT = m_ubxNavPVT;
    m_pos_ecef_0 = transformWGS84ToECEF(ubxNavPVT);
    m_R_ecefToLocal_0 = getR_ECEFToLocal(ubxNavPVT);
    printf("gps init values: %d, %f, %f, %f\r\n", GetNumSV(), m_pos_ecef_0(0), m_pos_ecef_0(1), m_pos_ecef_0(2));
}

NEOM9N::ubxNavPVT_t NEOM9N::GetUbxNavPVT()
{
    return m_ubxNavPVT;
}

// Position in ECEF in m
Eigen::Vector3f NEOM9N::GetPosECEF()
{
    return m_pos_ecef.cast<float>();
}

// Position in East-North-Up in m
Eigen::Vector3f NEOM9N::GetPosENU()
{
    return m_pos_enu.cast<float>();
}

// Velocity in East-North-Up in m
Eigen::Vector3f NEOM9N::GetVelENU()
{
    Eigen::Vector3f vel_enu;
    vel_enu << static_cast<float>( m_ubxNavPVT.velE ) * 1e-3f,
               static_cast<float>( m_ubxNavPVT.velN ) * 1e-3f,
              -static_cast<float>( m_ubxNavPVT.velD ) * 1e-3f;
    return vel_enu;
}

// 20  -    -   GNSSﬁx Type
uint8_t NEOM9N::GetFixType()
{
    return m_ubxNavPVT.fixType;
}

// 23  -    -   Number of satellites used in Nav Solution
uint8_t NEOM9N::GetNumSV()
{
    return m_ubxNavPVT.numSV;
}

//  0  -   ms   GPS time of week of the navigation epoch
uint8_t NEOM9N::GetGPSTimeOfWeek()
{
    return m_ubxNavPVT.iTOW;
}

// Heading of motion (2-D) in rad wrapped to (-pi, pi)
float NEOM9N::GetHeadMot()
{
    float headMot = static_cast<float>( m_ubxNavPVT.headMot )  * 1e-5f * M_PI / 180.0f;
    return atan2f( sinf(headMot), cosf(headMot) );
}

// 72 Heading accuracy estimate in rad (both motion and vehicle)
float NEOM9N::GetHeadAcc()
{
    return static_cast<float>( m_ubxNavPVT.headAcc )  * 1e-5f * M_PI / 180.0f;
}

// 40  Horizontal accuracy estimate in m
float NEOM9N::GethAcc()
{
    return static_cast<float>( m_ubxNavPVT.hAcc )  * 1e-3f;
}

// 44  Vertical accuracy estimate in m
float NEOM9N::GetvAcc()
{
    return static_cast<float>( m_ubxNavPVT.vAcc )  * 1e-3f;
}

// 68 Speed accuracy estimate in m/s
float NEOM9N::GetsAcc()
{
    return static_cast<float>( m_ubxNavPVT.sAcc )  * 1e-3f;
}

void NEOM9N::update()
{
    static char buffer[256];                // buffer to readout from the SerialBuffer
    static char check_buffer[4];            // buffer to check if a seuenz is starting
    const static int msg_buffer_len = 100;  // message length 94 + 2 (UBX_PAYLOAD_INDEX), do not know why +2, do not care
    static char msg_buffer[msg_buffer_len]; // buffer to decode message
    // ? | ? | PAYLOAD 92 | CK_A | CK_B
    static uint8_t msg_buffer_cntr = 0;
    static bool is_msg_sequenz = false;

    while(true) {

        ThisThread::flags_wait_any(threadFlag);

        if (m_bufferedSerial.readable()) {

            uint32_t read_len = m_bufferedSerial.read(buffer, sizeof(buffer));
            for (uint32_t i = 0; i < read_len; i++) {
                if (is_msg_sequenz) {
                    msg_buffer[msg_buffer_cntr++] = buffer[i];
                    if (msg_buffer_cntr == msg_buffer_len) { 
                        // checksum, we need to set those values because B5 62 01 07 are at the end of the message and we also need 01 and 07
                        uint8_t CK_A = 0x08; // 0x01 + 0x07
                        uint8_t CK_B = 0x09; // 0x01 + 0x01 + 0x07;
                        for (int i = 0; i < msg_buffer_len - 2 - 4; i++) {
                            CK_A += msg_buffer[i];
                            CK_B += CK_A;
                        }
                        // proceed if checksum is valid
                        if (CK_A == msg_buffer[msg_buffer_len - 2 - 4] && CK_B == msg_buffer[msg_buffer_len - 1 - 4]) {
                            ubxNavPVT_t ubxNavPVT = decodeUbxNavPVTmsg(msg_buffer);
                            // update fixType and numSV always
                            m_ubxNavPVT.fixType = ubxNavPVT.fixType;
                            m_ubxNavPVT.numSV = ubxNavPVT.numSV;
                            // update if fixtype 3 and numSV sufficent
                            if(ubxNavPVT.fixType == 3 && ubxNavPVT.numSV >= M_MIN_SATS) {
                                m_ubxNavPVT = ubxNavPVT;
                                static bool is_reset_local_internal = false;
                                if (!is_reset_local_internal) {
                                    reset_local();
                                    is_reset_local_internal = true;
                                }
                                m_pos_ecef = transformWGS84ToECEF(m_ubxNavPVT);
                                m_pos_enu = m_R_ecefToLocal_0 * ( m_pos_ecef - m_pos_ecef_0 );
                            }
                            is_msg_sequenz = false;
                            // HACK: update float sens_GPS[13];     // pos_enu(3), vel_enu(3), headMot, numSV, hAcc, vAcc, sAcc, headAcc, timestamp
                            /*
                            data.sens_GPS[0] = m_pos_enu(0);
                            data.sens_GPS[1] = m_pos_enu(1);
                            data.sens_GPS[2] = m_pos_enu(2);
                            Eigen::Vector3f vel_enu = GetVelENU();
                            data.sens_GPS[3] = vel_enu(0);
                            data.sens_GPS[4] = vel_enu(1);
                            data.sens_GPS[5] = vel_enu(2);
                            data.sens_GPS[6] = GetHeadMot();
                            data.sens_GPS[7] = (float)GetNumSV();
                            data.sens_GPS[8] = GethAcc();
                            data.sens_GPS[9] = GetvAcc();
                            data.sens_GPS[10] = GetsAcc();
                            data.sens_GPS[11] = GetHeadAcc();
                            data.sens_GPS[12] = global_timer.read();
                            */
#if PRINT_FOR_DEBUG
                            printf("%d, %d, %d, %d\r\n", ubxNavPVT.numSV, ubxNavPVT.lon, ubxNavPVT.lat, ubxNavPVT.height);
#endif                  
                        }
                    }
                }
                // shift past readings and check if a message starts
                check_buffer[3] = check_buffer[2]; // UBX_PVT_HEADER_0: 0xB5
                check_buffer[2] = check_buffer[1]; // UBX_PVT_HEADER_1: 0x62
                check_buffer[1] = check_buffer[0]; // UBX_PVT_CLASS   : 0x01
                check_buffer[0] = buffer[i];       // UBX_PVT_ID      : 0x07
                if(check_buffer[3] == UBX_PVT_HEADER_0 && check_buffer[2] == UBX_PVT_HEADER_1 && check_buffer[1] == UBX_PVT_CLASS && check_buffer[0] == UBX_PVT_ID) {
#if PRINT_FOR_DEBUG
                    int run_time = std::chrono::duration_cast<std::chrono::microseconds>(m_run_timer.elapsed_time()).count();
                    m_run_timer.reset();
                    printf("Time since last received package: %d, MSG len: %d\r\n", run_time, msg_buffer_cntr);
#endif
                    is_msg_sequenz = true;
                    msg_buffer_cntr = 0;
                }
            }
        }
    }
}

Eigen::Vector3d NEOM9N::transformWGS84ToECEF(const ubxNavPVT_t& ubxNavPVT)
{
    // persistent within this function
    const static double a = 6378137.0;         // WGS | 84 Earth semimajor axis
    const static double e = 8.1819191 * 1e-2;  // eccentricity

    double lon = static_cast<double>( ubxNavPVT.lon ) * 1e-7 * M_PI / 180.0;
    double lat = static_cast<double>( ubxNavPVT.lat ) * 1e-7 * M_PI / 180.0;
    double h   = static_cast<double>( ubxNavPVT.height ) * 1e-3;

    double sin_lat = sin(lat);
    double cos_lat = cos(lat);
    double N = a / sqrt(1 - e * e * sin_lat * sin_lat);

    Eigen::Vector3d pos_ecef;
    pos_ecef << (h + N) * cos_lat * cos(lon),
                (h + N) * cos_lat * sin(lon),
                (h + (1.0 - e*e)*N) * sin_lat;
    return pos_ecef;
}

Eigen::Matrix3d NEOM9N::getR_ECEFToLocal(const ubxNavPVT_t& ubxNavPVT)
{
    double lon = static_cast<double>(ubxNavPVT.lon) * 1e-7 * M_PI / 180.0;
    double lat = static_cast<double>(ubxNavPVT.lat) * 1e-7 * M_PI / 180.0;

    double sin_lat = sin(lat);
    double cos_lat = cos(lat);
    double sin_lon = sin(lon);
    double cos_lon = cos(lon);

    Eigen::Matrix3d R_ECEFToLocal;
    R_ECEFToLocal <<         -sin_lon,          cos_lon,       0,
                     -cos_lon*sin_lat, -sin_lat*sin_lon, cos_lat,
                      cos_lat*cos_lon,  cos_lat*sin_lon, sin_lat;
    return R_ECEFToLocal;
}

NEOM9N::ubxNavPVT_t NEOM9N::decodeUbxNavPVTmsg(const char *buf)
{
    static ubxNavPVT_t ubxNavPVT;
    static uint8_t index;

    /* example from u-center: len = 6*16+4 = 100 = 4 + 94 + 2
    09:41:53  0000  B5 62 01 07 5C 00 C8 72 AE 16 E6 07 06 09 09 29  µb..\.Èr®.æ....)
              0010  35 F3 FF FF FF FF 00 84 D7 17 00 00 A4 00 00 00  5óÿÿÿÿ..×...¤...
              0020  00 00 00 00 00 00 00 00 00 00 98 BD FF FF FF FF  ...........½ÿÿÿÿ
              0030  FF FF 00 34 F8 DF 00 00 00 00 00 00 00 00 00 00  ÿÿ.4øß..........
              0040  00 00 00 00 00 00 00 00 00 00 58 3E 0F 00 80 A8  ..........X>...¨
              0050  12 01 0F 27 00 00 EE 13 4F 2F 00 00 00 00 00 00  ...'..î.O/......
              0060  00 00 B9 53
    */

    // uint32_t iTOW;    //  0  -    -   GPS time of week of the navigation epoch
    index = UBX_PAYLOAD_INDEX + 0;
    ubxNavPVT.iTOW = buf[index++];
    ubxNavPVT.iTOW |= (buf[index++] << 8);
    ubxNavPVT.iTOW |= (buf[index++] << 16);
    ubxNavPVT.iTOW |= (buf[index++] << 24);
    // uint16_t year;    //  4  -    -   Year (UTC)
    index = UBX_PAYLOAD_INDEX + 4;
    ubxNavPVT.year = buf[index++];
    ubxNavPVT.year |= (buf[index++] << 8);
    // uint8_t month;    //  6  -    -   Month, range 1..12 (UTC)
    index = UBX_PAYLOAD_INDEX + 6;
    ubxNavPVT.month = buf[index++];
    // uint8_t day;      //  7  -    -   Day of month, range 1..31 (UTC)
    index = UBX_PAYLOAD_INDEX + 7;
    ubxNavPVT.day = buf[index++];
    // uint8_t hour;     //  8  -    -   Hour of day, range 0..23 (UTC)
    index = UBX_PAYLOAD_INDEX + 8;
    ubxNavPVT.hour = buf[index++];
    // uint8_t min;      //  9  -    -   Minute of hour, range 0..59 (UTC)
    index = UBX_PAYLOAD_INDEX + 9;
    ubxNavPVT.min = buf[index++];
    // uint8_t sec;      // 10  -    -   Seconds of minute, range 0..60 (UTC)
    index = UBX_PAYLOAD_INDEX + 10;
    ubxNavPVT.sec = buf[index++];
    // uint8_t fixType;  // 20  -    -   GNSSﬁx Type
    index = UBX_PAYLOAD_INDEX + 20;
    ubxNavPVT.fixType = buf[index];
    // uint8_t numSV;    // 23  -    -   Number of satellites used in Nav Solution
    index = UBX_PAYLOAD_INDEX + 23;
    ubxNavPVT.numSV = buf[index];
    // int32_t lon;      // 24 1e-7 deg  Longitude
    index = UBX_PAYLOAD_INDEX + 24;
    ubxNavPVT.lon = buf[index++];
    ubxNavPVT.lon |= (buf[index++] << 8);
    ubxNavPVT.lon |= (buf[index++] << 16);
    ubxNavPVT.lon |= (buf[index++] << 24);
    // int32_t lat;      // 28 1e-7 deg  Latitude
    index = UBX_PAYLOAD_INDEX + 28;
    ubxNavPVT.lat = buf[index++];
    ubxNavPVT.lat |= (buf[index++] << 8);
    ubxNavPVT.lat |= (buf[index++] << 16);
    ubxNavPVT.lat |= (buf[index++] << 24);
    // int32_t height;   // 32  -   mm   Height above ellipsoid
    index = UBX_PAYLOAD_INDEX + 32;
    ubxNavPVT.height = buf[index++];
    ubxNavPVT.height |= (buf[index++] << 8);
    ubxNavPVT.height |= (buf[index++] << 16);
    ubxNavPVT.height |= (buf[index++] << 24);
    // uint32_t hAcc;    // 40  -   mm   Horizontal accuracy estimate
    index = UBX_PAYLOAD_INDEX + 40;
    ubxNavPVT.hAcc = buf[index++];
    ubxNavPVT.hAcc |= (buf[index++] << 8);
    ubxNavPVT.hAcc |= (buf[index++] << 16);
    ubxNavPVT.hAcc |= (buf[index++] << 24);
    // uint32_t vAcc;    // 44  -   mm   Vertical accuracy estimate
    index = UBX_PAYLOAD_INDEX + 44;
    ubxNavPVT.vAcc = buf[index++];
    ubxNavPVT.vAcc |= (buf[index++] << 8);
    ubxNavPVT.vAcc |= (buf[index++] << 16);
    ubxNavPVT.vAcc |= (buf[index++] << 24);
    // int32_t velN;     // 48  -   mm/s NED north velocity
    index = UBX_PAYLOAD_INDEX + 48;
    ubxNavPVT.velN = buf[index++];
    ubxNavPVT.velN |= (buf[index++] << 8);
    ubxNavPVT.velN |= (buf[index++] << 16);
    ubxNavPVT.velN |= (buf[index++] << 24);
    // int32_t velE;     // 52  -   mm/s NED east velocity
    index = UBX_PAYLOAD_INDEX + 52;
    ubxNavPVT.velE = buf[index++];
    ubxNavPVT.velE |= (buf[index++] << 8);
    ubxNavPVT.velE |= (buf[index++] << 16);
    ubxNavPVT.velE |= (buf[index++] << 24);
    // int32_t velD;     // 56  -   mm/s NED down velocity
    index = UBX_PAYLOAD_INDEX + 56;
    ubxNavPVT.velD = buf[index++];
    ubxNavPVT.velD |= (buf[index++] << 8);
    ubxNavPVT.velD |= (buf[index++] << 16);
    ubxNavPVT.velD |= (buf[index++] << 24);
    // int32_t gSpeed;   // 60  -   mm/s Ground Speed (2-D)
    index = UBX_PAYLOAD_INDEX + 60;
    ubxNavPVT.gSpeed = buf[index++];
    ubxNavPVT.gSpeed |= (buf[index++] << 8);
    ubxNavPVT.gSpeed |= (buf[index++] << 16);
    ubxNavPVT.gSpeed |= (buf[index++] << 24);
    // int32_t headMot;  // 64 1e-5 deg  Heading of motion (2-D)
    index = UBX_PAYLOAD_INDEX + 64;
    ubxNavPVT.headMot = buf[index++];
    ubxNavPVT.headMot |= (buf[index++] << 8);
    ubxNavPVT.headMot |= (buf[index++] << 16);
    ubxNavPVT.headMot |= (buf[index++] << 24);
    // uint32_t sAcc;    // 68 -   mm/s Speed accuracy estimate
    index = UBX_PAYLOAD_INDEX + 68;
    ubxNavPVT.sAcc = buf[index++];
    ubxNavPVT.sAcc |= (buf[index++] << 8);
    ubxNavPVT.sAcc |= (buf[index++] << 16);
    ubxNavPVT.sAcc |= (buf[index++] << 24);
    // uint32_t headAcc; // 72 1e-5 deg  Heading accuracy estimate (both motion and vehicle)
    index = UBX_PAYLOAD_INDEX + 72;
    ubxNavPVT.headAcc = buf[index++];
    ubxNavPVT.headAcc |= (buf[index++] << 8);
    ubxNavPVT.headAcc |= (buf[index++] << 16);
    ubxNavPVT.headAcc |= (buf[index++] << 24);
    /*
    // int16_t magDec;   // 88 1e-2 deg  Magnetic declination.
    index = UBX_PAYLOAD_INDEX + 88;
    ubxNavPVT.magDec = buf[index++];
    ubxNavPVT.magDec |= (buf[index++] << 8);
    // uint16_t magAcc;  // 90 1e-2 deg  Magnetic declination accuracy
    index = UBX_PAYLOAD_INDEX + 90;
    ubxNavPVT.magAcc = buf[index++];
    ubxNavPVT.magAcc |= (buf[index++] << 8);
    */
    return ubxNavPVT;
}

void NEOM9N::sendThreadFlag()
{
    thread.flags_set(threadFlag);
}