#ifndef PGVL_H
#define PGVL_H

//! \brief Log an error internal to pgvl
#define LOGE(msg) \
   do { \
      std::cerr << "ERROR in pgvl: " << __func__ << "(): " << msg << std::endl; \
   } while(0); \

#endif /*PGVL_H*/
