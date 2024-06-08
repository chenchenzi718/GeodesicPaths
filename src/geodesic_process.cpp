#include "geodesic_process.h"
#include <igl/edge_lengths.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/exact_geodesic.h>
#include <corecrt_math_defines.h> //M_PI�ڴ�ͷ�ļ���
#include <igl/dijkstra.h>
#include <igl/adjacency_list.h>

GeodesicSolver::GeodesicSolver()
{

}

GeodesicSolver::~GeodesicSolver()
{

}


void GeodesicSolver::find_shortest_path(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXi& P)
{
    // ·���ϵ���ʼ�����յ�
    int start = P(0);
    std::set<int> end = { P(1) };

    // �����ڽӱ�
    std::vector<std::vector<int>> VV;
    igl::adjacency_list(F, VV);

    // �洢���·�������ǰ������
    Eigen::VectorXd min_distance;
    Eigen::VectorXi previous;

    // �������·��
    int found_target = igl::dijkstra(V, VV, start, end, min_distance, previous);

    // ����·��
    std::vector<int> path;
    path.clear();
    int current = P(1);
    while (current != start)
    {
        path.push_back(current);
        current = previous(current);
    }
    path.push_back(start);
    std::reverse(path.begin(), path.end());

    // ��·���洢�� P ��
    P.resize(path.size());
    for (int i = 0; i < path.size(); ++i)
    {
        P(i) = path[i];
    }
}

void GeodesicSolver::geodesic_remesh(const Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
    const Eigen::VectorXi& PathIndex, Eigen::MatrixXi& process_path)
{
    if (PathIndex.rows() < 3) {
        return;
    }

    for (int i = 1; i < PathIndex.rows() - 1; i++) {
        //  ȡ��·���ϵ����������㣬ȥ�����Ӧ��������Ƭ
        int a = PathIndex(i - 1);
        int b = PathIndex(i);
        int c = PathIndex(i + 1);

        //  Ѱ�� inner side һ��
        Eigen::MatrixXi T;
        bool is_tb = find_inner_triangles(V, F, a, b, c, T);

        // �� inner side һ��ִ�� flip edges
        flip_inner_triangles(V, F, a, b, c, is_tb, T, process_path);
    }
}


bool GeodesicSolver::flip_inner_triangles(const Eigen::MatrixXd& V, Eigen::MatrixXi& F,
    const int s, const int m, const int t, const bool is_tb,
    Eigen::MatrixXi& T, Eigen::MatrixXi& process_path)
{
    // �����һ�������Σ�������
    if (T.rows() < 2) {
        return false;
    }

    int index = 0;
    const int compare_amount = T.rows() - 1;
    for (int i = 0; i < compare_amount; i++) { 
        
        // ȡ��ÿ�����ڽӵ������Σ���ʼ���� angle �Ƿ�С�� pi
        int a = T(index, 0);
        int b = T(index, 1);
        int c = T(index, 2);
        int x = T(index + 1, 0);
        int y = T(index + 1, 1);
        int z = T(index + 1, 2);


        if (is_tb) {
            double opposite_angle = get_angle(V, a, b, c) + get_angle(V, x, z, y);

            if (opposite_angle > M_PI) {
                index += 1;
                continue;
            }

            // �޸Ĵ�ʱ·���ϵĵ� index ���� triangle �ϵ���Ϣ
            T(index, 1) = y;
            removeRow(T, index + 1);

            // ��һ������֮����һС�ε����·��Ϊ c->y
            process_path.conservativeResize(process_path.rows() + 1, 2);
            process_path(process_path.rows() - 1, 0) = c;
            process_path(process_path.rows() - 1, 1) = y;

            flip_edge(a, y, b, c, F);
        }
        else {
            double opposite_angle = get_angle(V, a, c, b) + get_angle(V, x, y, z);
            
            if (opposite_angle > M_PI) {
                index += 1;
                continue;
            }

            T(index, 2) = z;
            removeRow(T, index + 1);

            process_path.conservativeResize(process_path.rows() + 1, 2);
            process_path(process_path.rows() - 1, 0) = b;
            process_path(process_path.rows() - 1, 1) = z;

            flip_edge(a, b, c, z, F);
        }
    }
    return true;
}


bool GeodesicSolver::find_inner_triangles(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    const int a, const int b, const int c, Eigen::MatrixXi& T)
{
    // Ѱ�Ұ��� vertex b �����������Σ��洢�� T_AB ��
    Eigen::MatrixXi T_AB;
    for (int j = 0; j < F.rows(); j++)
    {
        int x = F(j, 0);
        int y = F(j, 1);
        int z = F(j, 2);
        if (b == x)
        {
            T_AB.conservativeResize(T_AB.rows() + 1, 3);
            T_AB.row(T_AB.rows() - 1) << x, y, z;
        }
        else if (b == y)
        {
            T_AB.conservativeResize(T_AB.rows() + 1, 3);
            T_AB.row(T_AB.rows() - 1) << y, z, x;
        }
        else if (b == z)
        {
            T_AB.conservativeResize(T_AB.rows() + 1, 3);
            T_AB.row(T_AB.rows() - 1) << z, x, y;
        }
    }

    // �����е������η�Ϊ�������֣�һ���Ǵ� b,a ����Ѱ�ҵ� b,a,x��Ȼ��� b,x ����Ѱ����һ�� b,x,y, ����Ѱ�ҵ�
    // b,z,c Ϊֹ����һϵ�е������γ�Ϊ T_A
    // �� b,c,x ��������Ѱ�ҵ��������γ�Ϊ T_B
    Eigen::MatrixXi T_A;
    Eigen::MatrixXi T_B;
    int temp = a;
    bool section_A = true;
    for (int j = 0; j < T_AB.rows(); j++)
    {
        for (int i = 0; i < T_AB.rows(); i++)
        {
            if (temp == T_AB(i, 1))
            {
                if (section_A)
                {
                    T_A.conservativeResize(T_A.rows() + 1, 3);
                    T_A.row(T_A.rows() - 1) << T_AB(i, 0), T_AB(i, 1), T_AB(i, 2);
                    temp = T_AB(i, 2);
                    if (temp == c)
                    {
                        section_A = false;
                    }
                }
                else
                {
                    if (temp == a)
                    {
                        break;
                    }
                    T_B.conservativeResize(T_B.rows() + 1, 3);
                    T_B.row(T_B.rows() - 1) << T_AB(i, 0), T_AB(i, 1), T_AB(i, 2);
                    temp = T_AB(i, 2);
                }
                continue;
            }
        }
    }

    // ͨ���Ƚ� T_A,T_B �ĽǶȺ����ж���һ���� a-b-c �����ߵ� inner side
    double angle_A = 0;
    for (int i = 0; i < T_A.rows(); i++)
    {
        int x = T_A(i, 0);
        int y = T_A(i, 1);
        int z = T_A(i, 2);
        angle_A += get_angle(V, z, x, y);
    }
    double angle_B = 0;
    for (int i = 0; i < T_B.rows(); i++)
    {
        int x = T_B(i, 0);
        int y = T_B(i, 1);
        int z = T_B(i, 2);
        angle_B += get_angle(V, z, x, y);
    }

    if (angle_A < angle_B)
    {
        T = T_A;
        return false;
    }
    else
    {
        T = T_B;
        // �� T_B �洢˳����һ�µ�������Ϊ�����Ǵ� b,c,x ���ϼ�¼�� b,z,a ��
        for (int i = 0; i < T.rows() / 2; i++)
        {
            T.row(i).swap(T.row(T.rows() - 1 - i));
        }
        return true;
    }
}


void GeodesicSolver::removeRow(
    Eigen::MatrixXi& matrix,
    unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

    matrix.conservativeResize(numRows, numCols);
}

int GeodesicSolver::get_triangle_index(const Eigen::MatrixXi& F, const int a, const int b, const int c)
{
    int index = -1;
    for (int i = 0; i < F.rows(); i++)
    {
        int x = F(i, 0);
        int y = F(i, 1);
        int z = F(i, 2);
        if (a == x && b == y && c == z)
        {
            index = i;
            break;
        }
        else if (a == y && b == z && c == x)
        {
            index = i;
            break;
        }
        else if (a == z && b == x && c == y)
        {
            index = i;
            break;
        }
    }
    return index;
}


bool GeodesicSolver::flip_edge(const int a, const int b, const int c, const int d, Eigen::MatrixXi& F)
{
    // ���������������Ϊ abc, acd
    //    b               b
    //   / \            / | \
    // a  -  c   ->   a   |   c
    //  \   /           \ |  /
    //    d               d

    // �½������� abd, bcd to F
    Eigen::MatrixXi newF(2, 3);
    newF << a, b, d,
        b, c, d;
    F.conservativeResize(F.rows() + 2, 3);
    F.block(F.rows() - 2, 0, 2, 3) = newF;

    // remove abc and bad
    removeRow(F, get_triangle_index(F, a, c, d));
    removeRow(F, get_triangle_index(F, a, b, c));

    return true;
}


double GeodesicSolver::get_angle(const Eigen::MatrixXd& V, const int a, const int b, const int c)
{
    // ���� bc ������ ab ����֮��ļн� arccos
    Eigen::Vector3d ab = V.row(a) - V.row(b);
    Eigen::Vector3d bc = V.row(c) - V.row(b);
    double angle = acos(ab.dot(bc) / (ab.norm() * bc.norm()));
    return angle;
}