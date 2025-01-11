import pandas as pd
import struct
import math
from gps_time import GPSTime  # gps-timeライブラリをインポート

# エンコードスケールファクタ（仮定値）
acc_encode_factor = 2**20 / (9.8 * 24)   # 加速度のスケールファクタ
R2D = 180 / math.pi
gyro_encode_factor = R2D / 4000.0 * 2**20  # 角速度のスケールファクタ

# Load the CSV file into a pandas DataFrame
csv_file_path = 'imu.csv'  # Update this with the correct path if necessary
imu_df = pd.read_csv(csv_file_path)

# Strip whitespace from column names
imu_df.columns = imu_df.columns.str.strip()

def encode_20bit_signed(value):
    return value & 0xFFFFF

def append_20bit_signed(buffer, values):
    """
    複数の20ビット符号付き整数をバッファに連結して追加する関数。
    buffer: バイト配列
    values: 符号付き20ビット整数のリスト
    """
    # valuesの長さが偶数であることを確認
    assert len(values) % 2 == 0, "Expected an even number of values for pairing"

    # 20ビットの符号付き整数をエンコードしてバッファに追加
    for i in range(0, len(values), 2):
        # ペアで20ビット符号付き整数を処理する
        encoded_value_0 = encode_20bit_signed(values[i])
        encoded_value_1 = encode_20bit_signed(values[i + 1])

        # 上位16ビットと下位4ビットに分割（1つ目の値）
        byte1_0 = (encoded_value_0 >> 12) & 0xFF
        byte2_0 = (encoded_value_0 >> 4) & 0xFF
        byte3_4bits_0 = encoded_value_0 & 0xF

        # 上位4ビットと下位16ビットに分割（2つ目の値）
        byte1_4bits_1 = (encoded_value_1 >> 16) & 0xF
        byte2_1 = (encoded_value_1 >> 8) & 0xFF
        byte3_1 = encoded_value_1 & 0xFF

        # 4ビットを結合
        combined_byte = (byte3_4bits_0 << 4) | byte1_4bits_1

        # バッファに5バイト分追加
        buffer.extend([byte1_0, byte2_0, combined_byte, byte2_1, byte3_1])

# Prepare the output list
imu_packet_list = []

for _, row in imu_df.iterrows():
    # GPS時刻（GPS TOW）を取得し、UNIX時刻に変換
    gps_tow = row['GPS TOW (s)']  # GPS TOWを秒として取得
    gps_week = row['GPS Week']     # GPS週を取得

    # GPS TOWとGPS Weekを使ってUNIX時刻に変換
    gps_time = GPSTime(week_number=gps_week, time_of_week=gps_tow).to_datetime()
    unix_time = gps_time.timestamp() + 32400
    """
    print("gps_week:{}", gps_week)
    print("gps_tow:{}", gps_tow)
    print("gps_time:{}", gps_time)
    print("unix_time:{}", unix_time)
    """

    # UNIX時刻から秒とナノ秒に分割
    sec = int(unix_time)  # 32400秒（9時間）を加算
    nsec = int((unix_time - sec) * 1e9)
    # その他のデータ処理
    acc_x = int(row['Acc X (m/s^2)'] * acc_encode_factor)
    acc_y = int(row['Acc Y (m/s^2)'] * acc_encode_factor)
    acc_z = int(row['Acc Z (m/s^2)'] * acc_encode_factor)
    gyro_x = int(row['Ang Rate X (deg/s)'] * gyro_encode_factor)
    gyro_y = int(row['Ang Rate Y (deg/s)'] * gyro_encode_factor)
    gyro_z = int(row['Ang Rate Z (deg/s)'] * gyro_encode_factor)

    # バッファの初期化
    preamble = 0xFECB
    reserved = 0x00
    buffer = bytearray()
    buffer.extend(preamble.to_bytes(2, 'big'))  # 16-bit preamble
    buffer.extend(reserved.to_bytes(1, 'big'))  # 8-bit reserved
    buffer.extend(struct.pack('>I', sec))       # 32-bit sec
    buffer.extend(struct.pack('>I', nsec))      # 32-bit nsec

    append_20bit_signed(buffer, [acc_x, acc_y, acc_z, gyro_x, gyro_y, gyro_z])

    data_length = len(buffer) - 3  # データ長を計算
    buffer[2] = data_length  # データ長をバッファに設定

    # CRCの計算と追加
    crc = 0
    for byte in buffer:
        crc ^= byte
    buffer.append(crc)

    # IMUパケットのリストに追加
    imu_packet_list.append(buffer)

# バイナリファイルとして保存
binary_file_path = 'imu.bin'
with open(binary_file_path, 'wb') as binary_file:
    for packet in imu_packet_list:
        binary_file.write(packet)

print(f'IMU packets saved to {binary_file_path}')