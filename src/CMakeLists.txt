SET( PGVL_SRCS
   BitmapImage.cpp
)

ADD_LIBRARY( pgvl
   SHARED
   ${PGVL_SRCS}
)