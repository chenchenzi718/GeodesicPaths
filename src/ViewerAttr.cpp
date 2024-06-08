#include "ViewerAttr.h"
#include<functional>
#include"geodesic_process.h"


ViewerAttr::ViewerAttr(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
    :V(V),F(F)
{
    viewer.append_mesh();
    F_his.push_back(F);

    set_mouse_down_action();
    set_key_press_action();
}


void ViewerAttr::set_mouse_down_action()
{
    viewer.callback_mouse_down =
        [&](igl::opengl::glfw::Viewer&, int, int) -> bool
    {
        if (vcs.CV.rows() >= 2) return false;
      
        // 从窗口坐标系转化到 opengl 坐标系
        last_mouse = Eigen::RowVector3f(
            viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
        if (vcs.placing_handles)
        {
            //  寻找鼠标点击的最近的 mesh 上的点
            int fid;
            Eigen::Vector3f bary;
            if (igl::unproject_onto_mesh(
                last_mouse.head(2),
                viewer.core().view,
                viewer.core().proj,
                viewer.core().viewport,
                V, F,
                fid, bary))
            {
                long c;
                bary.maxCoeff(&c);
                Eigen::RowVector3d new_c = V.row(F(fid, c));
               
                //  寻找找到的点是否在 vcs 内，如果不在就放入
                if (vcs.CV.rows() == 0 || (vcs.CV.rowwise() - new_c).rowwise().norm().minCoeff() > 0)
                {
                    // 将控制顶点放到 cv 内
                    vcs.CV.conservativeResize(vcs.CV.rows() + 1, 3);
                    vcs.CV.row(vcs.CV.rows() - 1) = new_c;
                    const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
                    draw_vertex(blue);

                    // 此时控制顶点 id 存放在 PathIndex 内
                    PathIndex.conservativeResize(PathIndex.rows() + 1);
                    PathIndex(PathIndex.rows() - 1) = F(fid, c);

                    return true;
                }
            }
        }
        return false;
    };
}


void ViewerAttr::draw_vertex(const Eigen::RowVector3d& color)
{
    viewer.data().set_points(vcs.CV, color);
}


void ViewerAttr::draw_edge(const Eigen::VectorXi& PathIndex, const Eigen::RowVector3d& color, 
    std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>>& edges)
{
    // 使用多次绘制的方式来加粗边
    Eigen::MatrixXd P1(PathIndex.size() - 1, 3);
    Eigen::MatrixXd P2(PathIndex.size() - 1, 3);
    for (int i = 0; i < PathIndex.size() - 1; ++i)
    {
        P1.row(i) = V.row(PathIndex(i));
        P2.row(i) = V.row(PathIndex(i + 1));
    }

    // 绘制原始边
    viewer.data().add_edges(P1, P2, Eigen::RowVector3d(color));
    // 记录绘制边
    edges.emplace_back(P1, P2, color);


    // 绘制加粗边
    for (double offset = 1e-6; offset <= 6e-6; offset += 1e-6) {
        Eigen::MatrixXd offset_P1 = P1;
        Eigen::MatrixXd offset_P2 = P2;
        offset_P1.col(0).array() += offset;
        offset_P2.col(0).array() += offset;
        viewer.data().add_edges(offset_P1, offset_P2, Eigen::RowVector3d(color));
        edges.emplace_back(offset_P1, offset_P2, color);

        offset_P1 = P1;
        offset_P2 = P2;
        offset_P1.col(0).array() -= offset;
        offset_P2.col(0).array() -= offset;
        viewer.data().add_edges(offset_P1, offset_P2, Eigen::RowVector3d(color));
        edges.emplace_back(offset_P1, offset_P2, color);

        offset_P1 = P1;
        offset_P2 = P2;
        offset_P1.col(1).array() += offset;
        offset_P2.col(1).array() += offset;
        viewer.data().add_edges(offset_P1, offset_P2, Eigen::RowVector3d(color));
        edges.emplace_back(offset_P1, offset_P2, color);

        offset_P1 = P1;
        offset_P2 = P2;
        offset_P1.col(1).array() -= offset;
        offset_P2.col(1).array() -= offset;
        viewer.data().add_edges(offset_P1, offset_P2, Eigen::RowVector3d(color));
        edges.emplace_back(offset_P1, offset_P2, color);
    }
}


void ViewerAttr::draw_edge(const std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>>& edges)
{
    // 是对上面的 draw_edge 方法的补足
    for (const auto& edge : edges) {
        viewer.data().add_edges(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge));
    }
}


void ViewerAttr::clear_drawn_point()
{
    viewer.data().clear_points();
}


void ViewerAttr::clear_drawn_edge()
{
    viewer.data().clear_edges();
}



void ViewerAttr::single_iteration()
{
    GeodesicSolver gsolver;

    if (vcs.CV.rows() > 1 && vcs.placing_handles) {
        // 寻找最短路径
        gsolver.find_shortest_path(V, F, PathIndex);

        // 连接选择出的控制点以及沿路的最短路径边
        const Eigen::RowVector3d red(1.0, 0., 0.);
        draw_edge(PathIndex, red, initial_shortest_edges);

        vcs.placing_handles = false;
        return;
    }

    //  执行两个点之间的 geodesic path 变换
    if (!vcs.placing_handles) {
        // 记录变化后的最短路径
        Eigen::MatrixXi process_path;

        if(!is_last_flip_change)
            gsolver.geodesic_remesh(V, F, PathIndex, process_path);

        if (F_his.back() != F)
        {
            // 说明依然在执行 geodesic 算法
            F_his.push_back(F);
            viewer.data().set_mesh(V, F);
        }
        else
        {
            is_last_flip_change = true;

            // 已经不用继续执行算法了，那么直接画出此时两个点之间的 geodesic path
            // 获取第一项和最后一项
            int first_element = PathIndex(0);
            int last_element = PathIndex(PathIndex.size() - 1);

            // 创建新的 PathIndex 并赋值
            Eigen::VectorXi newPathIndex(2);
            newPathIndex << first_element, last_element;
            gsolver.find_shortest_path(V, F, newPathIndex);

            const Eigen::RowVector3d blue(0., 1., 0.);
            clear_drawn_edge();
            draw_edge(newPathIndex, blue, geodesic_edges);
        }
    }
}


void ViewerAttr::roll_back()
{
    // 大于 1 才允许回退
    if (F_his.size() > 1)
    {
        F_his.pop_back();
        F = F_his.back();
        viewer.data().set_mesh(V, F);

        if (is_last_flip_change)
        {
            // 如果是最后一步回退，就将最后画的 geodesic 线先去掉，但是初始的最短路径不想去掉
            is_last_flip_change = false;
            clear_drawn_edge();
            draw_edge(initial_shortest_edges);
        }
    }
}


void ViewerAttr::clear_all_operation()
{
    F = F_his[0];
    F_his.clear();
    F_his.push_back(F);
    PathIndex.resize(0);

    initial_shortest_edges.clear();
    geodesic_edges.clear();
    is_last_flip_change = false;

    vcs.CV.resize(0, 3);  // 清空 CV 矩阵
    vcs.placing_handles = true;

    clear_drawn_point();
    clear_drawn_edge();

    viewer.data().set_mesh(V, F);
}


void ViewerAttr::set_key_press_action()
{
    viewer.callback_key_pressed =
        [&](igl::opengl::glfw::Viewer&, unsigned char key, int) -> bool
    {
        switch (key)
        {
        case ' ':
            single_iteration();
            break;
        case 'r':
            roll_back();
            break;
        case 'c':
            clear_all_operation();
            break;
        default:
            return false;
        }
        return true;
    };
}