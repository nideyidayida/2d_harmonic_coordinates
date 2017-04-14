//
//  Lasso.h
//  ex5
//
//  Created by Olga Diamanti on 21/04/15.
//
//

#ifndef __ex5__Lasso__
#define __ex5__Lasso__

#include <igl/embree/EmbreeIntersector.h>

//forward declaration of ViewerCore (needed for unprojection)
namespace igl {
  class ViewerCore;
}

class Lasso
{
public:
  
public:
  Lasso(const Eigen::MatrixXd &V_,
         const Eigen::MatrixXi &F_,
         const igl::ViewerCore &v);
  ~Lasso();
  void reinit();
private:
  const Eigen::MatrixXd &V;
  const Eigen::MatrixXi &F;
  igl::EmbreeIntersector ei;
  const igl::ViewerCore &viewercore;
  
  std::vector<std::vector<unsigned int> > stroke2DPoints;
  double d = -1;
public:
  int strokeAdd(int mouse_x, int mouse_y);
  void strokeFinish(Eigen::VectorXi &selected_vertices);
  int pickVertex(int mouse_x, int mouse_y);
  //the stroke
  std::vector< Eigen::Matrix<double, 1,3>  > strokePoints;
};

#endif /* defined(__ex5__Lasso__) */
