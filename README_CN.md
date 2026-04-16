# Picard 曲线枚举项目说明（中文版）

## 这个项目是做什么的？

枚举**小判别式的 Picard 曲线**。

Picard 曲线是一类亏格为 2 的代数曲线，形如：

```
y³ = f(x)
```

其中 `f(x)` 是四次多项式。这类曲线在数论中有特殊地位，判别式越小，曲线的算术结构越特殊（类似于椭圆曲线中的小判别式曲线）。

本项目有两个版本：
- **原版**：在有理数域 **Q** 上枚举，系数为整数
- **新版**：在 **Q(ζ₃)** 上枚举，系数为 Eisenstein 整数（形如 `a + b·ω`，其中 `ω = e^(2πi/3)`）

---

## 文件结构

```
mxm/
├── README.md              英文说明
├── README_CN.md           中文说明（本文件）
├── 1806.06289v5.pdf       相关论文
├── curves/                原版（Q 上枚举）
│   ├── pdisc_box.c        主程序
│   ├── mpoly128.c/h       多变量多项式树（加速判别式求值）
│   ├── mpolydisc128.c     判别式多项式的具体编码
│   ├── polyenum128.h      有限差分法（内层循环加速）
│   ├── mpzpolyutil.c      GMP 任意精度运算
│   ├── polyparse.c        多项式字符串解析
│   ├── cstd.c/h           工具函数
│   └── pdisc_box          编译好的可执行文件
├── curves_zeta3/          新版（Q(ζ₃) 上枚举）
│   ├── pdisc_box_zeta3.c  主程序（单文件）
│   ├── norm_of_disc.ipynb Sage 推导判别式公式的笔记本
│   └── pdisc_box_zeta3    编译好的可执行文件
└── results*.txt           运行结果
```

---

## 环境准备

所有操作都在 **WSL（Ubuntu）** 里进行，WSL 已安装在 `F:\WSL\Ubuntu\`。

打开 WSL 终端：按 `Win+R`，输入 `wsl`，回车。

### 依赖库

| 库 | 用途 | 安装命令 |
|---|---|---|
| GCC | C 编译器（含 OpenMP 支持） | 已预装 |
| libgmp-dev | 任意精度整数（原版需要） | `sudo apt-get install -y libgmp-dev` |

---

## 编译

### 原版（Q 上枚举）

```bash
cd ~/mxm/curves
gcc -O2 -fopenmp pdisc_box.c mpolydisc128.c mpoly128.c \
    mpzpolyutil.c polyparse.c cstd.c -lgmp -lm -o pdisc_box
```

### 新版（Q(ζ₃) 上枚举）

```bash
cd ~/mxm/curves_zeta3
gcc -O2 -fopenmp pdisc_box_zeta3.c -o pdisc_box_zeta3 -lm
```

> 两个可执行文件已编译好，直接跳到下面运行即可。

---

## 运行方法

### 原版：`pdisc_box`

```bash
~/mxm/curves/pdisc_box  系数上界  判别式上界  [线程数]
```

| 参数 | 含义 | 示例 |
|---|---|---|
| 系数上界 `c` | 搜索 \|fᵢ\| ≤ c 的所有系数组合 | `3` |
| 判别式上界 | 只输出 Δ ≤ 此值的曲线 | `1000000` |
| 线程数（可选） | 默认用全部 CPU 核心 | `4` |

**示例命令：**

```bash
# 系数 ≤ 3，判别式 ≤ 1,000,000，结果存入文件
~/mxm/curves/pdisc_box 3 1000000 > ~/mxm/results.txt

# 查看结果前 10 行
head -10 ~/mxm/results.txt

# 统计找到多少条曲线
wc -l ~/mxm/results.txt
```

**输出格式：**

```
判别式Δ:[多项式f(x)]
```

例如：
```
18915363:[x^4 - 3*x^2 - 3*x - 1]
```
表示曲线 `y³ = x⁴ - 3x² - 3x - 1`，其 Picard 曲线判别式 Δ = 18915363。

其中 Δ = 3⁹ · f₄³ · disc(f)²，disc(f) 是 f(x) 的多项式判别式（= 31²×3⁹ = 18915363）。

---

### 新版：`pdisc_box_zeta3`

```bash
~/mxm/curves_zeta3/pdisc_box_zeta3  系数上界  范数上界  [线程数]
```

| 参数 | 含义 | 示例 |
|---|---|---|
| 系数上界 `c` | 搜索 \|aᵢ\|, \|bᵢ\| ≤ c 的所有组合 | `3` |
| 范数上界 | 只输出 N(disc) ≤ 此值的曲线 | `1000000` |
| 线程数（可选） | 默认用全部 CPU 核心 | `4` |

**示例命令：**

```bash
# 系数 ≤ 3，判别式范数 ≤ 1,000,000，结果存入文件
~/mxm/curves_zeta3/pdisc_box_zeta3 3 1000000 > ~/mxm/results_zeta3.txt

# 查看结果前 10 行
head -10 ~/mxm/results_zeta3.txt

# 统计找到多少条曲线
wc -l ~/mxm/results_zeta3.txt
```

**输出格式：**

```
N(disc):[a0,b0,a1,b1,a2,b2]
```

例如：
```
18915363:[0,0,-1,0,0,0]
```

表示曲线 `y³ = x⁴ + (-1 + 0·ω)x + (0 + 0·ω)` = `y³ = x⁴ - x`，其判别式的范数为 18915363。

六个参数对应关系：

| 参数 | 含义 | 在曲线中的位置 |
|---|---|---|
| a0, b0 | 常数项 c₀ = a0 + b0·ω | `y³ = x⁴ + ... + c₀` |
| a1, b1 | 一次项系数 c₁ = a1 + b1·ω | `y³ = x⁴ + ... + c₁x + ...` |
| a2, b2 | 二次项系数 c₂ = a2 + b2·ω | `y³ = x⁴ + c₂x² + ...` |

完整曲线为：`y³ = x⁴ + (a2+b2·ω)x² + (a1+b1·ω)x + (a0+b0·ω)`

---

## 结果解读

### 原版结果（`results.txt`）

- 每行一条曲线，按找到的顺序排列（非排序）
- 判别式越小 → 曲线越特殊
- 可以用以下命令按判别式从小到大排序：

```bash
sort -t: -k1 -n ~/mxm/results.txt | head -20
```

### 新版结果（`results_zeta3.txt`）

- 每行一条 Eisenstein 整数系数的曲线
- `N(disc)` 是判别式从 Q(ζ₃) 到 Q 的范数，越小越特殊
- 排序查看：

```bash
sort -t: -k1 -n ~/mxm/results_zeta3.txt | head -20
```

---

## 速度参考

系数上界 `c = 3`，使用 8 线程：

| 版本 | 搜索规模 | 运行时间 | 找到曲线数 |
|---|---|---|---|
| 原版（Q） | `(2c+1)⁵ ≈ 1.7万` 条 | ~0.003 秒 | 约 700 条 |
| 新版（Q(ζ₃)） | `(2c+1)⁶ ≈ 11.8万` 条 | ~0.01 秒 | 约 2000 条 |

系数上界越大，搜索时间按 `c⁵`（原版）或 `c⁶`（新版）增长。

---

## 常用命令速查

```bash
# 进入 WSL
wsl

# 运行原版，系数 ≤ 5，判别式 ≤ 10⁷
~/mxm/curves/pdisc_box 5 10000000 > ~/mxm/results_c5.txt

# 运行新版，系数 ≤ 3，判别式范数 ≤ 10⁷
~/mxm/curves_zeta3/pdisc_box_zeta3 3 10000000 > ~/mxm/results_zeta3_c3.txt

# 按判别式排序，查看最小的 20 条
sort -t: -k1 -n ~/mxm/results_zeta3_c3.txt | head -20

# 查看运行日志（stderr）
~/mxm/curves_zeta3/pdisc_box_zeta3 3 10000000 > out.txt 2> log.txt
cat log.txt
```
