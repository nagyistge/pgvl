#ifndef POINT_H
#define POINT_H

/*!
 * \brief 2D integer coordinate
 */
class Point {
public:
   int x;
   int y;

   Point(int x = -1, int y = -1) :
      x(x),
      y(y)
   {
   }
   
   Point(Point const& other) :
      x(other.x),
      y(other.y)
   {
   }
   
   Point& operator=(Point const& rhs) {
      x = rhs.x;
      y = rhs.y;
      
      return *this;
   }
   
   bool operator==(Point const& rhs) const {
      return x == rhs.x && y == rhs.y;
   }
};

#endif /*POINT_H*/