#ifndef GEODESIC_PROCESS
#define GEODESIC_PROCESS
#include <Eigen/Core>


class GeodesicSolver {
public:
    GeodesicSolver();
    ~GeodesicSolver();

    // 利用 dijkstra 算法寻找最短路径，此时 P 内总共有两个点的 id
    void find_shortest_path(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXi& P);

    // 利用 flip_edge 算法得到两个点之间的 geodesic path
    void geodesic_remesh(const Eigen::MatrixXd& V, Eigen::MatrixXi& F, const Eigen::VectorXi& PathIndex,
        Eigen::MatrixXi& process_path);

#pragma region helper_function
    static void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);
    static int get_triangle_index(const Eigen::MatrixXi& F, const int a, const int b, const int c);
    static bool flip_edge(const int a, const int b, const int c, const int d, Eigen::MatrixXi& F);
    double get_angle(const Eigen::MatrixXd& V, const int a, const int b, const int c);
#pragma endregion

protected:
    bool flip_inner_triangles(const Eigen::MatrixXd& V, Eigen::MatrixXi& F,
        const int s, const int m, const int t, const bool is_tb,
        Eigen::MatrixXi& T, Eigen::MatrixXi& process_path);
    bool find_inner_triangles(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        const int a, const int b, const int c, Eigen::MatrixXi& T);
};

#endif