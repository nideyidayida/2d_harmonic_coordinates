#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/slice.h>
#include <igl/triangle/triangulate.h>
#include <igl/project.h>
#include <igl/unproject.h>


//activate this for alternate UI (easier to debug)
//#define UPDATE_ONLY_ON_UP

using namespace std;

//vertex array, #V x3
Eigen::MatrixXd V3D(0,3);
//vertex array, #V x2
Eigen::MatrixXd V(0,2);
//face array, #F x3
Eigen::MatrixXi F(0,3);

//cage vertex array, #C x2
Eigen::MatrixXd C(0,2);
//cage edge arary, #C x2
Eigen::MatrixXi E(0,2);
//vertex array after triangulation, #V2 x2
Eigen::MatrixXd V2;
//face array after triangulation, #F2 x3
Eigen::MatrixXi F2;
//cage function, #V2-#C x#C
Eigen::MatrixXd Hfunction(0,0);

//mouse interaction
enum MouseMode { SELECT, TRANSLATE, NONE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;

//index of cage vertex being moved
int moving_cage = -1;

//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

//function declarations (see below for implementation)
bool solve(igl::viewer::Viewer& viewer);
Eigen::MatrixXd computeH();
void addCageVertex(Eigen::RowVector2d cageP);
int getMovnCage(Eigen::RowVector2d cageP);
Eigen::RowVector2d pickPoint(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y);
Eigen::MatrixXd readMatrix(const char *filename);

bool callback_mouse_down(igl::viewer::Viewer& viewer, int button, int modifier);
bool callback_mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(igl::viewer::Viewer& viewer, int button, int modifier);
bool callback_pre_draw(igl::viewer::Viewer& viewer);
bool callback_key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers);

bool solve(igl::viewer::Viewer& viewer)
{
  Eigen::MatrixXd newInternalP = Hfunction * C;
  for (int v = 0; v < V.rows(); v++) {
    V.row(v) = newInternalP.row(v);
  }
  
  V3D.col(0) = V.col(0); V3D.col(1) = V.col(1);
  
  return true;
};

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << "Usage project mesh.off" << endl;
    exit(0);
  }

  // Read mesh
  igl::readOFF(argv[1],V3D,F);
  assert(V3D.rows() > 0);

  V.resize(V3D.rows(), 2);
  V.col(0) = V3D.col(0); V.col(1) = V3D.col(1);

  // Plot the mesh
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = callback_key_down;
  viewer.callback_init = [&](igl::viewer::Viewer& viewer)
  {

      viewer.ngui->addGroup("Deformation Controls");

      viewer.ngui->addVariable<MouseMode >("MouseMode",mouse_mode)->setItems({"SELECT", "TRANSLATE", "NONE"});

//      viewer.ngui->addButton("ClearSelection",[](){ selected_v.resize(0,1); });

      viewer.screen->performLayout();
      return false;
  };

  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  viewer.callback_pre_draw = callback_pre_draw;

  viewer.data.clear();
  viewer.data.set_mesh(V3D, F);

  viewer.core.point_size = 10;
  viewer.core.set_rotation_type(igl::viewer::ViewerCore::ROTATION_TYPE_TRACKBALL);

  viewer.launch();
}

void addCageVertex(Eigen::RowVector2d cageP) {
  C.conservativeResize(C.rows()+1, 2);
  C.row(C.rows()-1) = cageP;
  
  // connect cage vertices.
  E.resize(C.rows(), 2);
  for (int c = 0; c < C.rows(); c++) {
    E(c, 0) = c;
    E(c, 1) = (c+1) % C.rows();
  }
}

Eigen::MatrixXd computeH() {
  Eigen::MatrixXd CV(C.rows()+V.rows(), 2);
  CV << C, V;
  Eigen::MatrixXd H;
  igl::triangle::triangulate(CV, E, H, "a100.0.05q", V2, F2);
  // compute bi-Laplacian coefficient.
  Eigen::SparseMatrix<double> L, M, Minv;
  igl::cotmatrix(V2, F2, L);
  igl::massmatrix(V2, F2, igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::invert_diag(M, Minv);
  Eigen::SparseMatrix<double> A = L * Minv * L;

  // extract Aff and Afc from A.
  Eigen::MatrixXi internalP(V2.rows() - C.rows(), 1);
  Eigen::MatrixXi cageC(C.rows(), 1);
  for (int v = 0; v < V2.rows()-C.rows(); v++) {
    internalP(v, 0) = v + C.rows();
  }
  for (int c = 0; c < C.rows(); c++) {
    cageC(c, 0) = c;
  }
  Eigen::SparseMatrix<double> Aff, Afc;
  igl::slice(A, internalP, internalP, Aff);
  igl::slice(A, internalP, cageC, Afc);

  // solve for H.
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(Aff);
  return solver.solve(-1*Eigen::MatrixXd(Afc));
}

int getMovnCage(Eigen::RowVector2d cageP) {
  int movingCage = 0;
  double distance = (C.row(0) - cageP).norm();
  for (int c = 1; c < C.rows(); c++) {
    double disTmp = (C.row(c) - cageP).norm();
    if (disTmp < distance) {
      distance = disTmp;
      movingCage = c;
    }
  }
  return movingCage;
}

Eigen::RowVector2d pickPoint(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y) {
  // Cast a ray in the view direction starting from the mouse position
  double x = mouse_x;
  double y = viewer.core.viewport(3) - mouse_y;  
  Eigen::RowVector3d pt;  
  Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
  // project selected point onto viewer window.
  Eigen::Vector3f proj = igl::project(pt.transpose().cast<float>().eval(), modelview, viewer.core.proj,viewer.core.viewport);
  float c = proj[2];
  // unproject the point based on the object.
  pt = igl::unproject(Eigen::Vector3f(x, y, c), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();
  
  Eigen::RowVector2d pt2d;
  pt2d(0) = pt(0); pt2d(1) = pt(1);

  return pt2d;
}

bool callback_mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
  if (button == (int) igl::viewer::Viewer::MouseButton::Right)
    return false;

  down_mouse_x = viewer.current_mouse_x;
  down_mouse_y = viewer.current_mouse_y;

  if (mouse_mode == SELECT) {
    addCageVertex(pickPoint(viewer, viewer.current_mouse_x, viewer.current_mouse_y));
    doit = true;
  } else if (mouse_mode == TRANSLATE) {
    moving_cage = getMovnCage(pickPoint(viewer, viewer.current_mouse_x, viewer.current_mouse_y));
    doit = true;
  }
  return doit;
}

bool callback_mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y)
{
  if (!doit)
    return false;

  if (mouse_mode == TRANSLATE) {
    C.row(moving_cage) = (pickPoint(viewer, viewer.current_mouse_x, viewer.current_mouse_y));

#ifndef UPDATE_ONLY_ON_UP
    solve(viewer);
    down_mouse_x = mouse_x;
    down_mouse_y = mouse_y;
#endif
    return true;
  }
  return false;
}

bool callback_mouse_up(igl::viewer::Viewer& viewer, int button, int modifier)
{
  if (!doit)
    return false;
  doit = false;

  if (mouse_mode == TRANSLATE) {
#ifdef UPDATE_ONLY_ON_UP
    if(moving_cage>=0)
      solve(viewer);
#endif
    moving_cage = -1;

    return true;
  }
  return false;
};

bool callback_pre_draw(igl::viewer::Viewer& viewer)
{
  // Initialize vertex colors
  vertex_colors = Eigen::MatrixXd::Constant(V.rows(),3,.9);

  //clear points and lines
  viewer.data.set_points(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
  viewer.data.set_edges(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXi::Zero(0,3), Eigen::MatrixXd::Zero(0,3));

  // update the vertex position all the time
  viewer.data.V.resize(V3D.rows(),3);
  viewer.data.V << V3D;

  viewer.data.dirty |= viewer.data.DIRTY_POSITION;

  Eigen::MatrixXd C3d(C.rows(), 3);
  C3d.col(0) = C.col(0); C3d.col(1) = C.col(1);
  viewer.data.add_points(C3d, Eigen::MatrixXd::Constant(C.rows(),3,.2));
  viewer.data.set_edges(C3d, E, Eigen::MatrixXd::Constant(C.rows(),3,.2));

  return false;
}

bool callback_key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
  bool handled = false;
  if (key == 'H') {
    Hfunction = computeH();
    callback_key_down(viewer, '1', 0);
    handled = true;
  }
  if (key == 'S')
  {
    mouse_mode = SELECT;
    callback_key_down(viewer, '1', 0);
    handled = true;
  }

  if (key == 'N')
  {
    mouse_mode = NONE;
    callback_key_down(viewer, '1', 0);
    handled = true;
  }

  if (key == 'M') {
    mouse_mode = TRANSLATE;
    callback_key_down(viewer, '1', 0);
    handled = true;
  }

  viewer.ngui->refresh();
  return handled;
}