
import pandas as pd
import struct
import math

# エンコードスケールファクタ（仮定値）
acc_encode_factor = 2**20 / (9.8 * 24)   # 加速度のスケールファクタ
R2D = 180 / math.pi
gyro_encode_factor = R2D / 4000.0 * 2**20  # 角速度のスケールファクタ
# Load the CSV file into a pandas DataFrame
csv_file_path = 'sample.csv'  # Update this with the correct path if necessary
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
    pending_bits = None  # 前回の4ビットを保持する変数

    encoded_value_0 = encode_20bit_signed(values[0])
    # 上位16ビットと下位4ビットに分割
    byte1_0 = (encoded_value_0 >> 12) & 0xFF   # 上位8ビット
    byte2_0 = (encoded_value_0 >> 4) & 0xFF    # 中間8ビット
    byte3_4bits_0 = encoded_value_0 & 0xF      # 下位4ビット

    encoded_value_1 = encode_20bit_signed(values[1])
    # 上位4ビットと下位16ビットに分割
    byte1_4bits_1 = (encoded_value_1 >> 16) & 0xF  # 上位4ビット
    byte2_1 = (encoded_value_1 >> 8) & 0xFF        # 中間8ビット
    byte3_1 = encoded_value_1 & 0xFF               # 下位8ビット

    # 下位4ビットと上位4ビットを結合
    combined_byte = (byte3_4bits_0 << 4) | byte1_4bits_1

    # バッファに追加
    buffer.extend([byte1_0, byte2_0, combined_byte, byte2_1, byte3_1])

    encoded_value_2 = encode_20bit_signed(values[2])
    # 上位16ビットと下位4ビットに分割
    byte1_2 = (encoded_value_2 >> 12) & 0xFF   # 上位8ビット
    byte2_2 = (encoded_value_2 >> 4) & 0xFF    # 中間8ビット
    byte3_4bits_2 = encoded_value_2 & 0xF      # 下位4ビット

    encoded_value_3 = encode_20bit_signed(values[3])
    # 上位4ビットと下位16ビットに分割
    byte1_4bits_3 = (encoded_value_3 >> 16) & 0xF  # 上位4ビット
    byte2_3 = (encoded_value_3 >> 8) & 0xFF        # 中間8ビット
    byte3_3 = encoded_value_3 & 0xFF               # 下位8ビット

    # 下位4ビットと上位4ビットを結合
    combined_byte_2 = (byte3_4bits_2 << 4) | byte1_4bits_3

    # バッファに追加
    buffer.extend([byte1_2, byte2_2, combined_byte_2, byte2_3, byte3_3])

    encoded_value_4 = encode_20bit_signed(values[4])
    # 上位16ビットと下位4ビットに分割
    byte1_4 = (encoded_value_4 >> 12) & 0xFF   # 上位8ビット
    byte2_4 = (encoded_value_4 >> 4) & 0xFF    # 中間8ビット
    byte3_4bits_4 = encoded_value_4 & 0xF      # 下位4ビット

    encoded_value_5 = encode_20bit_signed(values[5])
    # 上位4ビットと下位16ビットに分割
    byte1_4bits_5 = (encoded_value_5 >> 16) & 0xF  # 上位4ビット
    byte2_5 = (encoded_value_5 >> 8) & 0xFF        # 中間8ビット
    byte3_5 = encoded_value_5 & 0xFF               # 下位8ビット

    # 下位4ビットと上位4ビットを結合
    combined_byte_3 = (byte3_4bits_4 << 4) | byte1_4bits_5

    # バッファに追加
    buffer.extend([byte1_4, byte2_4, combined_byte_3, byte2_5, byte3_5])

# Prepare the output list
imu_packet_list = []

for _, row in imu_df.iterrows():
    # Extract values from the DataFrame row
    sec = int(row['GPS TOW (s)'])  # assuming the GPS TOW is used for seconds (32 bits unsigned)
    nsec = int((row['GPS TOW (s)'] - sec) * 1e9)  # deriving nanoseconds part (32 bits unsigned)

    # Accelerations and angular velocities are 20-bit signed integers

    acc_x = int(row['Acc X (m/s^2)'] * acc_encode_factor)
    acc_y = int(row['Acc Y (m/s^2)'] * acc_encode_factor)
    acc_z = int(row['Acc Z (m/s^2)'] * acc_encode_factor)
    gyro_x = int(row['Ang Rate X (deg/s)'] * gyro_encode_factor)
    gyro_y = int(row['Ang Rate Y (deg/s)'] * gyro_encode_factor)
    gyro_z = int(row['Ang Rate Z (deg/s)'] * gyro_encode_factor)

    # Initialize a buffer with preamble and reserved bytes (16-bit preamble + 8-bit reserved)
    preamble = 0xFECB
    reserved = 0x00

    # Create a byte buffer with preamble and reserved field
    buffer = bytearray()
    buffer.extend(preamble.to_bytes(2, 'big'))  # 16-bit preamble
    buffer.extend(reserved.to_bytes(1, 'big'))  # 8-bit reserved

    # Append the data fields (32-bit unsigned for sec and nsec)
    buffer.extend(struct.pack('>I', sec))       # 32-bit sec
    buffer.extend(struct.pack('>I', nsec))      # 32-bit nsec

    append_20bit_signed(buffer, [acc_x, acc_y, acc_z, gyro_x, gyro_y, gyro_z])

    # Set message length without header and parity
    # Length is calculated as the length of the data fields (excluding preamble and reserved)
    data_length = len(buffer) - 3  # 3 bytes are preamble and reserved
    buffer[2] = data_length  # Set length in the reserved field position

    # Calculate CRC-8 of the entire message (excluding the CRC byte itself)
    crc = 0
    for byte in buffer:
        crc ^= byte

    # Append CRC to the message
    buffer.append(crc)

    # Append the final buffer to the list of IMU packets
    imu_packet_list.append(buffer)

# Save the final generated IMU packets in binary format
binary_file_path = 'sample.bin'

with open(binary_file_path, 'wb') as binary_file:
    for packet in imu_packet_list:
        binary_file.write(packet)

print(f'IMU packets saved to {binary_file_path}')
