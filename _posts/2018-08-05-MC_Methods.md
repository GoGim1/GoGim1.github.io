---
layout:     post
title:      "蒙特卡罗方法的数学基础"
subtitle:   "Mathematical Foundations of Monte Carlo Methods"
date:       2018-08-05 01:00:00
author:     "Wu"
header-img: "img/temp.jpg"
tags:
    - 计算机图形学
    - 翻译
    - 读书笔记
---

>本文是我在 [scratchapixel](www.scratchapixel.com) 的学习笔记


## 简介

### 快速介绍

生活中的许多事情都难以准确评估，特别是涉及到非常大的数字的时候。例如计算一个国家的成年人口的平均身高，我们不太可能去统计每一个人的身高，而是取一个人口样本并计算其平均身高，得出其结果作为近似值。虽然结果并不十分精确，但随着样本的大小增加，近似值与实际结果之间的误差随着样本量的增加而变小。这是一种用准确性换取时间的方法。
![Approximation(Average(X)) = { 1 \over N} \sum_{n=1}^N x_n.](http://latex.codecogs.com/gif.latex?Approximation(Average(X))&space;=&space;{&space;1&space;\over&space;N}&space;\sum_{n=1}^N&space;x_n.)
注意，这里的X是一个**随机变量**（如前面每个人的高度），在统计中，随机变量X的平均值称为**期望值**，写为E(X)。
![](http://latex.codecogs.com/gif.latex?E(X)&space;\approx&space;{&space;1&space;\over&space;N&space;}&space;\sum_{n&space;=&space;1}^N&space;x_n.)

### 运用蒙特卡罗方法的光线追踪

在计算机中，图片的每一个像素都被保存为单一颜色，而这个颜色应该由每个像素「看到」的全部颜色共同决定的。如图，![图1](http://www.scratchapixel.com/images/upload/monte-carlo-methods/areacam1.png?)
理想情况下，每个像素的颜色应该这么计算：![2](http://latex.codecogs.com/gif.latex?L_{pixel}&space;=&space;\int_{pixel&space;area}&space;L(x_p)&space;dA,)

然而，我们可以通过在像素区域中选择几个随机样本位置并平均它们的颜色，求得每个像素的颜色:![3](http://www.scratchapixel.com/images/upload/monte-carlo-methods/areacam2.png?)
公式如下：![4](http://latex.codecogs.com/gif.latex?L_{pixel}&space;\approx&space;{1&space;\over&space;N&space;}&space;\sum_{n=1}^N&space;L(x_n),)

这种方法是有缺陷的。首先它的结果并不是精确的，其次是每次计算都会得到不同的结果。我们可以通过增加样本数量来克服这些问题。但是，要将误差减少一半，你需要两倍的样本。换句话说，它们的收敛速度（随着样本数量的增加，它们收敛到正确结果的速度）非常低。

![](http://www.scratchapixel.com/images/upload/monte-carlo-methods/mcintegration01.png?)
上面是随机取8个样本的结果。

![](http://www.scratchapixel.com/images/upload/monte-carlo-methods/mcintegration03.png?)
上面是随机取样，每种进行16次计算的结果。

每次计算结果的差异程度也取决于网格颜色之间的变化程度。如图：![](http://www.scratchapixel.com/images/upload/monte-carlo-methods/mcintegration04.png?)我们可以看到其计算的差异程度远远大于前一个例子。

我们用于计算的代码如下：
```
#include <fstream> 
#include <cstdlib> 
#include <cstdio> 
 
int main(int argc, char **argv) 
{ 
        std::ifstream ifs; 
        ifs.open("./tex.pbm"); 
        std::string header; 
        uint32_t w, h, l; 
        ifs >> header; 
        ifs >> w >> h >> l; 
        ifs.ignore(); 
        unsigned char *pixels = new unsigned char[w * h * 3]; 
        ifs.read((char*)pixels, w * h * 3); 
        // sample
        int nsamples = 8; 
        srand48(13); 
        float avgr = 0, avgg = 0, avgb = 0; 
        float sumr = 0, sumg = 0, sumb = 0; 
        for (int n = 0; n < nsamples; ++ n) { 
                float x = drand48() * w; 
                float y = drand48() * h; 
                int i = ((int)(y) * w + (int)(x)) * 3; 
                sumr += pixels[i]; 
                sumg += pixels[i + 1]; 
                sumb += pixels[i + 2]; 
        } 
        sumr /= nsamples; 
        sumg /= nsamples; 
        sumb /= nsamples; 
        for (int y = 0; y < h; ++y) { 
                for (int x = 0; x < w; ++x) { 
                        int i = (y * w + x) * 3; 
                        avgr += pixels[i]; 
                        avgg += pixels[i + 1]; 
                        avgb += pixels[i + 2]; 
                } 
        } 
        avgr /= w * h; 
        avgg /= w * h; 
        avgb /= w * h; 
        printf("Average %0.2f %0.2f %0.2f, Approximation %0.2f %0.2f %0.2f\n", avgr, avgg, avgb, sumr, sumg, sumb); 
        delete [] pixels; 
        return 0; 
} 
```

如果随着样本大小增加，其结果收敛到预期的值，则估计值是**无偏**的。结果收敛到其他的值，则估计值是**有偏**的。看起来，无偏估计似乎比有偏估计更好。但如果偏差足够小并且有偏估计比无偏估计收敛得更快，又或者有偏估计比无偏估计的方差小，那么你可以考虑选择有偏估计。

![](http://www.scratchapixel.com/images/upload/monte-carlo-methods/mcintegration05.png?)
假如我们要得到从摄像机看到的P点的颜色，我们需要计算照射到P点上半球面的光线。![](http://latex.codecogs.com/gif.latex?L_P&space;=&space;\int_\Omega&space;L_\Omega.)

同样的，我们可以选择在半球上取样，以获得到达P的光的近似值。![](http://latex.codecogs.com/gif.latex?L_P&space;\approx&space;{&space;1&space;\over&space;N&space;}&space;\sum_{n=1}^N&space;L_n.)

这个过程是可以递归的,我们使用递归过程来计算到达可见点P的光。![](http://www.scratchapixel.com/images/upload/monte-carlo-methods/mcintegration06.png?)
代码如下：
```
Vec3f monteCarloIntegration(P, N) 
{ 
    Vec3f lightAtP = 0; 
    int nsamples = 8; 
    for (int n = 0; n < nsamples; ++n) { 
        Ray sampleRay = sampleRayAboveHemisphere(P, N); 
        Vec3f Phit; 
        Vec3f Nhit; 
        if (traceRay(sampleRay, Phit, Nhit)) { 
            lightAtP += monteCarloIntegration(Phit, Nhit); 
        } 
    } 
    lightAtP /= nsamples; 
 
     return lightAtP; 
} 
 
void render() 
{ 
    Ray r; 
    computeCameraRayDir(r, ...); 
    Vec3f Phit; 
    Vec3f Nhit; 
    if (traceRay(r, Phit, Nhit)) { 
        Vec3f lightAtP = monteCarloIntegration(Phit, Nhit); 
        ... 
    } 
} 
```

## 随机变量和概率

## 概率分布：第1部分
## 概率属性
## 统计学概论
## 期望值
## 方差和标准差
## 概率分布：第2部分
## 采样分布
## 概率密度函数（PDF）和累积分布函数（CDF）
## 随机变量函数的期望值：无意识统计学家的规律
## 逆变换采样方法
## 估计



