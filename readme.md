# 实现 geodesic paths 的框架说明

简要说明一下框架用法。

## External Libraries

* [Eigen](https://gitlab.com/libeigen/eigen/-/releases/3.4-rc1), version 3.4.0.
* [libigl](https://github.com/libigl/libigl).
* [GLAD](https://glad.dav1d.de/), version 4.6.
* [GLFW](https://www.glfw.org/), version 3.3.9.

但是只需要使用的时候另外下载 Eigen 即可，其他库都已经下载并放在文件夹 /3rd_party 内。需要修改 CmakeLists 内的 Eigen 路径为你的下载路径，在文件的第 31 行。
```
SET(EIGEN_PATH "D:/myapp/Eigen/eigen-3.4.0")
```


## Usage

修改完 CmakeLists 之后，直接用 Cmake GUI 进行 cmake generate，build 得到 C++ 项目。打开 GeodesicPaths.sln 解决方案，将工程 GeodesicPaths 设置为启动项目点击运行即可。

得到的 GUI 界面可以用鼠标左键选择 mesh 上的两个点，然后点击 `空格键` 寻找到这两个点之间的最短路径，以红色标注。然后`连续多次`点击 `空格键` 就可以不断的进行 geodesic paths 的算法，最终将两个点之间的 geodesic paths 以绿色进行标注。

如果想要回退到上一次计算，点击 `r` 键即可；如果想要清除所有的选择结果，点击 `c` 键即可。
