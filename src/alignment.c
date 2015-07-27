
#include "alignment.h"

/* main function. */
int main(int argc, char *argv[]) {
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	printf("%s\n", argv[argc-1]);
	ks1->s = "ACTGTACAACTGAAGGACTaaaGACATGGCAATCCTTAAGAATTTTACCTACAGAATGAATGCACACATATATCATTCCAAATTTGTCAACTTATATAAATATGGTCCCATTTTCACAGTTAATTGGCTCAGTACAGTCACCCATAAGCCACTGTTCTTTAAACAGGAAGCTAAATTTATTTACATAGGAAGCTGCAATTTATTCCCATCTCAATATGGCAATAAATGGAGATATACAAGCTTTGTAAAGAAAAAGAGATCCATACGTACTGGGGAAATACTCATGTGTTTCATGCAAAAGAATTATATTCATAAGGAAACCA";
	ks2->s = "ATCCACTTTTACAATATGTAAAAGGTACTTTTAACTTCCTTTCATTGAACCAGTGTACAACAGTTCACTGTACAACTGAAGGACTGACATGGCAATCCTTAAGAATTTTACCTACAGAATGAATGCACACAatttgacacattttcttagtttcaaaagattatttaaaaaaggaattcagtagattgacttgtaaataaccattgcagattttgaatctgcaaaaatccgtcacattgctgttgggacagattaagataaggctaaaatttttttccaagttattgcaaactgctatggaaagagaaagtaatccaaaatgtataaaatggcccatggacaaatccaaaccacgcaatttttgtaaataaaggtttattgcaatatggccacatctacttactcatgTATATCATTCCAAATTTGTCAACTTATATAAATATGGTCCCATTTTCACAGTTAATTGGCTTCACCAAGTAAGAAAATATGGGTAAAAACACAATTCAAGGTCACTCAAGTTTATCATCCTCGTAAGTAACAACAGCTCTCTATTTGAAGGTATATGGGAATCTCAAGTAGAATATTCAAGACTTTCTTAACAATATGTAAATAGCTCAATAATAGCAAACTTTAAAATGTACATTTTTCTCAGTACAGTCACCCATAAGCCACTGTTCTTTAAACAGGAAGCTAAATTTATTTACATAGGAAGCTGCAATTTATTCCCATCTCAATATGGCAATAAATGGAGATATACAAGTAttcgctatattaggtaagttatatccgataggacaggtttacatttacagttactcaattcatgttaataagcttttttccagCCATGTTTCAACTACTTTGTAAAGAAAAAGAGATCCATACGTACTGGGGAAATACTCATGTGTTTCATGTATTTTTCAAAAGAATTATATTCATAAGGAAACCACCACTTAATTATTCCAATCTCATGCAGGCTCAATGATTTTTTAAGTCAGTATTTTGTCCTATGGAACTAATCCAACCGAGTAACAGATGATTTTAAAAATACGCTCATTTAAATACAAACA";
	ks1->l = strlen(ks1->s);
	ks2->l = strlen(ks2->s);
	//kstring_read(argv[argc-1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	if(ks1->l > ks2->l) die("first sequence must be shorter than the second\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align(ks1, ks2, r1, r2));
	printf("%s\n%s\n", r1->s, r2->s);
	//kstring_destory(ks1);
	//kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
