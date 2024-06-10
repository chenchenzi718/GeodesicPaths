#include "geodesic_process.h"
#include"ViewerAttr.h"

#include <igl/read_triangle_mesh.h>
#include <string>
#include <iostream>


int main(int argc, char* argv[])
{
    const Eigen::RowVector3d grey(0.5, 0.5, 0.5); // 模型设置为灰色

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh("../example/bunny_10k.obj", V, F);

    ViewerAttr viewer_attr = ViewerAttr(V, F);

    viewer_attr.viewer.data().set_mesh(viewer_attr.V, viewer_attr.F);
    viewer_attr.viewer.data().set_colors(grey);
    viewer_attr.viewer.data().show_vertex_labels = true;

    viewer_attr.viewer.launch();

    return EXIT_SUCCESS;
}