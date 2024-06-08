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
      
        // �Ӵ�������ϵת���� opengl ����ϵ
        last_mouse = Eigen::RowVector3f(
            viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
        if (vcs.placing_handles)
        {
            //  Ѱ�������������� mesh �ϵĵ�
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
               
                //  Ѱ���ҵ��ĵ��Ƿ��� vcs �ڣ�������ھͷ���
                if (vcs.CV.rows() == 0 || (vcs.CV.rowwise() - new_c).rowwise().norm().minCoeff() > 0)
                {
                    // �����ƶ���ŵ� cv ��
                    vcs.CV.conservativeResize(vcs.CV.rows() + 1, 3);
                    vcs.CV.row(vcs.CV.rows() - 1) = new_c;
                    const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
                    draw_vertex(blue);

                    // ��ʱ���ƶ��� id ����� PathIndex ��
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
    // ʹ�ö�λ��Ƶķ�ʽ���Ӵֱ�
    Eigen::MatrixXd P1(PathIndex.size() - 1, 3);
    Eigen::MatrixXd P2(PathIndex.size() - 1, 3);
    for (int i = 0; i < PathIndex.size() - 1; ++i)
    {
        P1.row(i) = V.row(PathIndex(i));
        P2.row(i) = V.row(PathIndex(i + 1));
    }

    // ����ԭʼ��
    viewer.data().add_edges(P1, P2, Eigen::RowVector3d(color));
    // ��¼���Ʊ�
    edges.emplace_back(P1, P2, color);


    // ���ƼӴֱ�
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
    // �Ƕ������ draw_edge �����Ĳ���
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
        // Ѱ�����·��
        gsolver.find_shortest_path(V, F, PathIndex);

        // ����ѡ����Ŀ��Ƶ��Լ���·�����·����
        const Eigen::RowVector3d red(1.0, 0., 0.);
        draw_edge(PathIndex, red, initial_shortest_edges);

        vcs.placing_handles = false;
        return;
    }

    //  ִ��������֮��� geodesic path �任
    if (!vcs.placing_handles) {
        // ��¼�仯������·��
        Eigen::MatrixXi process_path;

        if(!is_last_flip_change)
            gsolver.geodesic_remesh(V, F, PathIndex, process_path);

        if (F_his.back() != F)
        {
            // ˵����Ȼ��ִ�� geodesic �㷨
            F_his.push_back(F);
            viewer.data().set_mesh(V, F);
        }
        else
        {
            is_last_flip_change = true;

            // �Ѿ����ü���ִ���㷨�ˣ���ôֱ�ӻ�����ʱ������֮��� geodesic path
            // ��ȡ��һ������һ��
            int first_element = PathIndex(0);
            int last_element = PathIndex(PathIndex.size() - 1);

            // �����µ� PathIndex ����ֵ
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
    // ���� 1 ���������
    if (F_his.size() > 1)
    {
        F_his.pop_back();
        F = F_his.back();
        viewer.data().set_mesh(V, F);

        if (is_last_flip_change)
        {
            // ��������һ�����ˣ��ͽ���󻭵� geodesic ����ȥ�������ǳ�ʼ�����·������ȥ��
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

    vcs.CV.resize(0, 3);  // ��� CV ����
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