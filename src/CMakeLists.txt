SET( PGVL_SRCS
   BitmapImage.cpp
   ImageProcessing.cpp
   ppm.cpp
)

ADD_LIBRARY( pgvl
   SHARED
   ${PGVL_SRCS}
)