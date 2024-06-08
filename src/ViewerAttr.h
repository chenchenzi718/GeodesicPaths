#ifndef VIEWER_ATTR
#define VIEWER_ATTR

#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <Eigen/Core>

class ViewerAttr
{
public:
	ViewerAttr(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	~ViewerAttr() {};

	void set_mouse_down_action();			// ֻ����ѡ��������
	void set_key_press_action();            // ������������һ��ѭ��

	void draw_vertex(const Eigen::RowVector3d& color);
	void draw_edge(const Eigen::VectorXi& PathIndex, const Eigen::RowVector3d& color,
		std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>>& edges);
	void draw_edge(const std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>>& edges);
	void clear_drawn_point();
	void clear_drawn_edge();
	
	void single_iteration();             // ��Ӧ��һ��ѭ��
	void roll_back();                    // ��Ӧ��һ�λ��˲���
	void clear_all_operation();            // �����еĽ��ȫ������


	Eigen::MatrixXd V;              // ��� mesh �� vertex����СΪ #V*3
	Eigen::MatrixXi F;				// ��� mesh �� face��ÿ�� face ������ vertex id����˴�С�� #F*3 
	Eigen::VectorXi PathIndex;      // ���ѡ��� path �ϵ� vertex id
	std::vector<Eigen::MatrixXi> F_his;  // ��� flip �仯�� F ��ʷ

	// ��¼���ǻ滭���ʼ�����·���Լ����ջ滭�����·��
	std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>> initial_shortest_edges, geodesic_edges;
	bool is_last_flip_change = false;    // ��¼���һ�� flip edge �任������֮�󲻻����޸� F ����ʷ


	// ����ָʾ��ǰ viewer �ϵ�ѡ��״̬
	struct VertexChooseState
	{
		Eigen::MatrixXd CV;         // n*3 ��С�Ŀ��Ƶ�����ܹ�ѡ��������
		bool placing_handles = true;
	} vcs;

	Eigen::RowVector3f last_mouse;  // ��¼���ѡ������һ����
	igl::opengl::glfw::Viewer viewer;  
};

#endif