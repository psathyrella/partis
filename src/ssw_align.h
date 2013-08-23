#ifndef SSW_ALIGN_H
#define SSW_ALIGN_H

void align_reads (const char* ref_path, const char* qry_path,
                  const char* output_path,
                  const int32_t match, /* 2 */
                  const int32_t mismatch, /* 2 */
                  const int32_t gap_open, /* 3 */
                  const int32_t gap_extension); /* 1 */

#endif
