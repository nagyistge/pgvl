SET( PGVL_SRCS
   Image.cpp
   ImageProcessing.cpp
   ppm.cpp
)

ADD_LIBRARY( pgvl
   SHARED
   ${PGVL_SRCS}
)

TARGET_LINK_LIBRARIES( pgvl ${SDL2_LIBRARY} )
