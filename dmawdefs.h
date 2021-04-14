#include <sdsl/bit_vectors.hpp>
#include <sdsl/cst_sct3.hpp>
typedef sdsl::bp_interval<> node_type;

#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"
#define INPUT_STR               "input_str.txt"

#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYVOUBZJX*"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

using namespace sdsl;
using namespace std;

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

typedef tuple<INT,INT,INT,INT> mytuple;
struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_filename;
   unsigned int         total_length;
 };
 
 struct TMaw
 {
   INT	letter;
   INT	pos;
   INT 	size;
 };

double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );
unsigned int RevComStr ( unsigned char * str, unsigned char * str2, INT iLen );
unsigned int compute_dmaw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, unsigned char * seq1, unsigned char * seq1_id, unsigned int twoseq );
unsigned int compute_maw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, TMaw ** Occ, INT * NOcc );
unsigned char Mapping( int a );
int RevMapping ( unsigned char b );
unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );

unsigned int GetBefore (
				unsigned char * seq,
                                INT n,
                                int sigma,
				INT * SA,
                                INT * LCP,
                                bit_vector * Before,
                                bit_vector * Beforelcp );
unsigned int GetMaws(
				unsigned char * seq,
				unsigned char * seq_id,
				INT * SA,
				INT n,
				int sigma,
				INT * LCP,
				bit_vector * Before,
				bit_vector * Beforelcp,
				unsigned int k,
				unsigned int K,
				char * out_file,
				TMaw ** Occ, INT * NOcc );
