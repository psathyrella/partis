#ifndef IG_ALIGN_H
#define IG_ALIGN_H

#include <stdint.h>

void ig_align_reads(const char *ref_path,
                    const uint8_t n_extra_refs,
                    const char **extra_ref_paths,
                    const char *qry_path,
                    const char *output_path,
                    const int32_t match,     /* 2 */
                    const int32_t mismatch,  /* 2 */
                    const int32_t gap_o,     /* 3 */
                    const int32_t gap_e,     /* 1 */
                    const unsigned max_drop, /* 1000 */
                    const unsigned bandwidth,
                    const uint8_t n_threads,
                    const char *read_group,
                    const char *read_group_id);

#endif
