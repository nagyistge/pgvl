SET( PGVL_SRCS
   BitmapImage.cpp
   ImageProcessing.cpp
)

ADD_LIBRARY( pgvl
   SHARED
   ${PGVL_SRCS}
)