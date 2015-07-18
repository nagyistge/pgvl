#ifndef POINT_H
#define POINT_H

/*!
 * \brief 2D integer coordinate
 */
class Point {
public:
   //! \brief X (horizontal) coordinate
   int x;
   //! \brief Y (vertical) coordinate
   int y;

   //! \brief Default constructor
   Point(int x = -1, int y = -1) :
      x(x),
      y(y)
   {
   }

   //! \brief Copy constructor
   Point(Point const& other) :
      x(other.x),
      y(other.y)
   {
   }

   //! \brief Assignment operator
   Point& operator=(Point const& rhs) {
      x = rhs.x;
      y = rhs.y;
      
      return *this;
   }

   //! \brief Equality operator
   bool operator==(Point const& rhs) const {
      return x == rhs.x && y == rhs.y;
   }
};

#endif /*POINT_H*/