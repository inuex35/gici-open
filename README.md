# 高精度測位チャレンジ2024用 GICI-LIB README

## GNSS/INS/Camera 統合ナビゲーションライブラリ

GICI-LIBは、GNSS、INS（慣性航法システム）、カメラを統合したナビゲーションのためのオープンソースソフトウェアです。本ライブラリを用いて、高精度測位チャレンジ2024向けのデータ処理およびアルゴリズム実験が可能です。

以下の構成では、東京および名古屋のデータセットを使用した各種設定ファイルについて説明し、インストールから実行までの手順を詳述します。

---

## 1. 必要環境

### 1.1 OS
- Ubuntu 20.04（推奨）
- Ubuntu 18.04 または 22.04 でも動作確認済み

### 1.2 必須ライブラリ
1. **Eigen**
   - バージョン: 3.3 以上
   - 線形代数用のC++テンプレートライブラリ。
   - [Eigen公式ページ](https://eigen.tuxfamily.org/index.php?title=Main_Page)

2. **OpenCV**
   - バージョン: 4.2.0 以上
   - コンピュータビジョンライブラリ。
   - [OpenCV公式ページ](https://opencv.org/releases/)

3. **Yaml-cpp**
   - バージョン: 0.6.0 以上
   - YAML形式のデコード/エンコードライブラリ。
   - [Yaml-cpp公式ページ](https://github.com/jbeder/yaml-cpp)

4. **Glog と Gflags**
   - バージョン: 0.6.0 以上
   - ロギング制御ライブラリ。
   - [Glog公式ページ](https://github.com/google/glog)
   - [Gflags公式ページ](https://github.com/gflags/gflags)

5. **Ceres-Solver**
   - バージョン: 2.1.0 以上
   - 非線形最適化ライブラリ。
   - [Ceres-Solver公式ページ](http://ceres-solver.org/)

---

## 2. インストール手順

### 2.1 通常ビルド

```bash
cd <gici-root-directory>
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```

GICI-LIBを実行するには以下のコマンドを使用してください。

```bash
./build/gici_main <gici-config-file>
```

---

## 3. 設定ファイル

以下の設定ファイルを使用して各地域のデータを処理できます。

- **東京**
  - `option/tc.yaml` ：trainのtokyo1 の設定ファイル
  - `option/tc1.yaml`：testのtokyo1 の設定ファイル
  - `option/tc2.yaml`：testのtokyo2 の設定ファイル
  - `option/tc3.yaml`：testのtokyo3 の設定ファイル

- **名古屋**
  - `option/nagoya_tc.yaml` ：trainのnagoya1 の設定ファイル
  - `option/nagoya1_tc.yaml`：testのnagoya1 の設定ファイル
  - `option/nagoya2_tc.yaml`：testのnagoya2 の設定ファイル
  - `option/nagoya3_tc.yaml`：testのnagoya3 の設定ファイル

---

## 4. 実行手順

### 4.1 東京データセットの場合

1. 設定ファイルを指定してGICI-LIBを実行します。

   ```bash
   ./build/gici_main ./option/tc1.yaml
   ```

2. 必要に応じて他の設定ファイル（例：`tc2.yaml`, `tc3.yaml`）を使用してください。

### 4.2 名古屋データセットの場合

1. 設定ファイルを指定してGICI-LIBを実行します。

   ```bash
   ./build/gici_main ./option/nagoya1_tc.yaml
   ```

2. 必要に応じて他の設定ファイル（例：`nagoya2_tc.yaml`, `nagoya3_tc.yaml`）を使用してください。

---

## 5. データセットのダウンロード

高精度測位チャレンジ用データセットは以下からダウンロード可能です。

[[PPC-Dataset](https://github.com/taroz/PPC-Dataset)
]

---

## 6. 可視化

### 非ROS環境での可視化

以下はで[RTKLIB](https://rtklib.com/)を用いた可視化の例です。
RTKLIBで可視化ができるのでconfig内のポートに接続してください。

![無題の動画 ‐ Clipchampで作成 (9)](https://github.com/user-attachments/assets/4d102374-5d52-4d9a-945b-a22a24862eae)


---

## 7. ライセンス

GICI-LIBは[GPL v3](https://www.gnu.org/licenses/gpl-3.0.html)のもとで配布されています。改変および配布は自由ですが、GPL v3の条件を満たす必要があります。
