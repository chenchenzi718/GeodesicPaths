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

	void set_mouse_down_action();			// 只允许选择两个点
	void set_key_press_action();            // 设置怎样进行一次循环

	void draw_vertex(const Eigen::RowVector3d& color);
	void draw_edge(const Eigen::VectorXi& PathIndex, const Eigen::RowVector3d& color,
		std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>>& edges);
	void draw_edge(const std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>>& edges);
	void clear_drawn_point();
	void clear_drawn_edge();
	
	void single_iteration();             // 对应着一次循环
	void roll_back();                    // 对应着一次回退操作
	void clear_all_operation();            // 将所有的结果全部回退


	Eigen::MatrixXd V;              // 存放 mesh 的 vertex，大小为 #V*3
	Eigen::MatrixXi F;				// 存放 mesh 的 face，每个 face 有三个 vertex id，因此大小是 #F*3 
	Eigen::VectorXi PathIndex;      // 存放选择的 path 上的 vertex id
	std::vector<Eigen::MatrixXi> F_his;  // 存放 flip 变化的 F 历史

	// 记录我们绘画的最开始的最短路径以及最终绘画的最短路径
	std::vector<std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::RowVector3d>> initial_shortest_edges, geodesic_edges;
	bool is_last_flip_change = false;    // 记录最后一次 flip edge 变换，在这之后不会再修改 F 的历史


	// 用来指示当前 viewer 上的选择状态
	struct VertexChooseState
	{
		Eigen::MatrixXd CV;         // n*3 大小的控制点矩阵，总共选择两个点
		bool placing_handles = true;
	} vcs;

	Eigen::RowVector3f last_mouse;  // 记录鼠标选择的最后一个点
	igl::opengl::glfw::Viewer viewer;  
};

#endif