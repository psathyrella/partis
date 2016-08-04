#ifndef SSW_ALIGN_H
#define SSW_ALIGN_H

#include <stdint.h>

void align_reads(const char *ref_path,
                 const char *qry_path,
                 const char *output_path,
                 const int32_t match,    /* 2 */
                 const int32_t mismatch, /* 2 */
                 const int32_t gap_o,    /* 3 */
                 const int32_t gap_e,    /* 1 */
                 const uint8_t n_threads,
                 const int32_t n_keep,
                 const int32_t max_drop,
                 const char *read_group,
                 const char *read_group_id);

#endif
