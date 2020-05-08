# 密度泛函理论计算液晶相变

## 编译

将`dft.cpp`编译两次
```shell
mpic++ -std=c++11 -O3 dft.cpp -o pre
mpic++ -std=c++11 -O3 dft.cpp -o main
```

## 运行
第一次用MC产生需要的能量表
```shell
mpirun -np 16 ./pre 10 30
```
表示计算长度为10，内转角为30的FCh分子


第二次使用迭代计算ODF，并计算性质
```shell
./main 10 30