#ifndef MEEP_MT_H
#define MEEP_MT_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void meep_mt_init_genrand(unsigned long s);
void meep_mt_restore_genrand();
void meep_mt_init_by_array(unsigned long init_key[], int key_length);
unsigned long meep_mt_genrand_int32(void);
long meep_mt_genrand_int31(void);
double meep_mt_genrand_real1(void);
double meep_mt_genrand_real2(void);
double meep_mt_genrand_real3(void);
double meep_mt_genrand_res53(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MEEP_MT_H */
