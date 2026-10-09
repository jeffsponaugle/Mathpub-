#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

// 
// A19 (1350168131352524443) occurs at 8,858,170,606 and 4,750,253,204 digits, 
// A20 (84756845106452435773) occurs at 17601613331 and 15494062638 digits.
// A21 (585270898631522188621) occurs at 69,658,148,238 - 5,644,388,308  (earliest)
//      685151818447609104295 = 685151818447609104295
//      Locations: 95,098,922,518 - 2,126,170,845
//      953657151565326368763 = 953657151565326368763
//      Locations: 99,532,319,752 - 61,194,089,794
// A22 (2761994111668451704865) occurs at 122,027,419,722 - 73,134,869,029
//      1587124362959858673273 = 1587124362959858673273 
//      Locations: 132,915,426,432 - 107,265,929,645
//      Defines for size/hashtable
//      2761994111668451704865 = 2761994111668451704865 (earliest)
//      Locations: 122,027,419,722 - 73,134,869,029

#define CANSIZE  24                      // Size of char buffer for candidate numbers
#define HASHSIZE 50000000000             // Size of hashtable (50,000,000,000) bytes, 400 billion bits
#define PISOURCE "../pi-10billionX.txt"

// A197123
// The core structure used for each candidate number.

// Global variables ( used for performance )
int8_t *S_TopBloom;                  //  Pointer to hash table
uint64_t S_TotalElements = 0;
uint64_t S_TotalElementsLast = 0;          //  Total number of elements added to array.
uint64_t S_UnallocatedHash = 0;
long S_digitsize;                       //  The size 'digit match' we are looking for
struct timeval S_stop, S_start, S_last;         //  Start and stop time metrics.
long S_MaxDepth = 0;
long S_MaxDepthCount = 0;
char *S_NewBlockPtr=0;
long S_RemainingElementsinBlock=0;
char *S_NextBlockPtr;
uint64_t S_BloomPass = 0;
uint64_t qtr[100] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                      2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                      3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3
                   };

FILE *datafile;
FILE *pisrc;
char S_digitbuffer[4194304];            // (1024 x 4096)
int S_digitbuffercount = 0;

void Complete();
void AddNumberToBloom(char* num,uint64_t offset);
long totaltimer();
long lasttimer();
void Cleanup();
char getdigitbuffer();


int main(int argc, char *argv[])
{
        
        long digitfilter = 0;           // If we are filtering, this is the starting digit.
        int filterenabled = 0;          // 1 == filtering enabled

        printf("\r\n\r\nStarting PI Search\r\n");
        // Read input parameters
        // seqhash <DIGITSTOTESTinM> <SIZEDIGITTOLOOKFOR> {FILTERNUMBER}
        char *p;
        // counter is the number of digits of PI we are going to search.
        int64_t counter = (uint64_t)strtol(argv[1], &p, 10) * (1000000);

        // S_digitsize is the number of digits we are searching for.
        S_digitsize = strtol(argv[2], &p, 10);

        if (S_digitsize < 19)
        {
                printf("Digit length must be 23 or more.\r\n");
                exit(0);
        }

        // If the filter arguement is present, that tells us to limit to numbers
        //  that start with the specified digit (/10 calculations)
        if (argc==4)
        {
           digitfilter = strtol(argv[3], &p, 10); 
           printf("Limiting to [%ld]\r\n",digitfilter);
           filterenabled=1;
        }

        char filename[128];
        if ( filterenabled == 1)
                sprintf(filename,"pidata_%ld_%lld_F%ld.txt",S_digitsize,counter,digitfilter);
        else
                sprintf(filename,"pidata_%ld_%lld.txt",S_digitsize,counter);

        printf("Datalog file name:%s\r\n",filename);
        datafile = fopen(filename,"w");
        printf("Processing %3.2f M numbers. \r\n",(float)counter/1000000.0);
        printf("Sequence Length %ld\r\n",S_digitsize);
        printf("Allocating bloom table\r\n");

        // Allocate the hashtable.
        S_TopBloom = malloc(((uint64_t)sizeof(int8_t))*(uint64_t)HASHSIZE);
        
        printf("Allocated bloom table of size: %lld\r\n",((uint64_t)sizeof(int8_t))*(uint64_t)HASHSIZE );
        printf("Buffer location: %lx\r\n",(long unsigned int)S_TopBloom);

        char num[CANSIZE];              // num is the number we are currently testing.

        printf("Clearing bloomtable space.\r\n");
        // Clear out all hash table entries
	memset(S_TopBloom,0,HASHSIZE);
        S_UnallocatedHash = (uint64_t)HASHSIZE*8;
        // clear the testing number.
        for (int u=0; u<CANSIZE;u++)
                num[u]=0;        

        printf("Starting calculations.\r\n");
        gettimeofday(&S_start, NULL);

        // ppi-10billion.txt contains the first 10 billion digits of pi, starting with
        //  3.14159265358.... 
        // We will skip over the '3.'
        printf("Opening PISOURCE\r\n");
        pisrc = fopen(PISOURCE,"rb");
        if (pisrc == 0)
        {
                printf("Error opening PISOURCE %s\r\n");
                exit(1);
        }
        getdigitbuffer();              // read the 3
        getdigitbuffer();              // read the .

        // Read in S_digitsize digits for the first number to test.
        for (int z=0; z<S_digitsize;z++)
        {
                num[z] = getdigitbuffer();
        }

        uint64_t tries = 0;
        char nextdigit;

        do 
        {       
                // Insert this number into the array, cheching first to see if it should be filtered.  
                if (filterenabled==1)
                {
                        if (num[0] == ((char)digitfilter + '0'))
                        {
                           S_TotalElements++; 
                           AddNumberToBloom(num,tries+1);  
                          
                        }
                }
                else
                {
                    S_TotalElements++;
                    AddNumberToBloom(num,tries+1);   
                 
                }

                // At this point the number has beed added, so if it makes it here it was not matched.
                // Every 1% of total tries we print a status report

                if (tries%(counter/100)==0)
                {
                        float hashratio = ((float)S_UnallocatedHash / ((float) HASHSIZE*8))*100.0;
			printf("\r\nCurrent Element Count: %6.2fM : %5.2fM  PassCount %lld (%2.5f) (Unallocated %lld(%2.4f%%))   - \r\n",S_TotalElements/1000000.0,(float)(tries/1000000.0),S_BloomPass,((float)S_BloomPass/(float)S_TotalElements),S_UnallocatedHash,hashratio);

                         long diff = totaltimer()+1;  // diff in ms
                         long difflast = lasttimer()+1;
                         float rate = S_TotalElements / (uint64_t)diff;  // rate in elem/ms, or kelem/s
                         float ratelast = (S_TotalElements-S_TotalElementsLast) / (uint64_t)difflast;
                         printf("       Cumulative Rate: %4.2f K Elements/sec   Section Rate: %4.2f K Elements/sec\r\n",rate,ratelast);
                         S_TotalElementsLast = S_TotalElements;
                         fflush(datafile);

                        //getrusage(RUSAGE_SELF,&r_usage);
                        // Print the maximum resident set size used (in kilobytes).
                        //printf("Memory usage: %6.2f MB\n",r_usage.ru_maxrss/1000000.0);
                }

                // get the next digit, then move all the digits over by 1.

                memcpy(&num[0],&num[1],(S_digitsize-1));
                //for (int t=1; t<S_digitsize; t++)
                //        num[t-1]=num[t];
                num[S_digitsize-1] = getdigitbuffer();
                tries++;
                
        } while (tries < counter);

        Complete();
        //getrusage(RUSAGE_SELF,&r_usage);
        // Print the maximum resident set size used (in kilobytes).
        //printf("Memory usage: %6.2f MB\n",r_usage.ru_maxrss/1000000.0);
        //PrintCandidateListSize();
        //PrintCandidateList();
        
}

char getdigitbuffer()
{
        if (S_digitbuffercount == 0)
        {
                fread(S_digitbuffer,4194304,1,pisrc);
                S_digitbuffercount = 4194304;
        }
        return S_digitbuffer[4194304-S_digitbuffercount--];
}


void Cleanup()
{
    uint64_t sum =0;

    fclose(datafile);
    printf("Calculating Bloom Table Statistics\r\n");
    {
        uint64_t id;
        for (id=0;id<HASHSIZE;id++)
        {
                sum +=  __builtin_popcount(S_TopBloom[id]);
        }
        printf("Total Bits Used: %llu (%2.6f)\r\n",sum, (double)sum/((double)HASHSIZE*8.0));
    }
}

void AddNumberToBloom(char* num, uint64_t offset)
{
        short int BloomFlag1,BloomFlag2,BloomFlag3, BloomFlag4;
        int r;

        // Lets start by figuring out which hash entry to use.    We used each digit in the 
        // number to pick a hash value.  The last value is used to select between the first 
        // and second 10B hash values. 
    
        uint64_t hashentry1 = (((uint64_t)(num[1] - '0'))) 
                        + (uint64_t)(num[2] - '0')  * 10 
                        + (uint64_t)(num[3] - '0')  * 100 
                        + (uint64_t)(num[4] - '0')  * 1000
                        + (uint64_t)(num[5] - '0')  * 10000
                        + (uint64_t)(num[6] - '0')  * 100000
                        + (uint64_t)(num[7] - '0')  * 1000000
                        + (uint64_t)(num[8] - '0')  * 10000000                // 100M
                        + (uint64_t)(num[9] - '0')  * 100000000 
                        + (uint64_t)(num[10] - '0') * 1000000000
                        + (uint64_t)(num[11] - '0') * 10000000000     // 10b
                        + ((qtr[(num[12] - '0')* 10 + (num[13] - '0')]) * (uint64_t)100000000000);
                                                                                              
        if ((hashentry1 > (uint64_t)3999999999999))
        {
                printf("Bad Hash1 Value: %llu, %s\r\n",hashentry1,num);
                exit(0);
        }
        
	uint64_t hashentry2 =  ((((uint64_t)(num[12] - '0')))) 
                        + (uint64_t)(num[13] - '0') * 10 
                        + (uint64_t)(num[14] - '0') * 100 
                        + (uint64_t)(num[15] - '0') * 1000
                        + (uint64_t)(num[16] - '0') * 10000
                        + (uint64_t)(num[17] - '0') * 100000
                        + (uint64_t)(num[18] - '0') * 1000000
                        + (uint64_t)(num[9] - '0')  * 10000000                // 100M
                        + (uint64_t)(num[0] - '0')  * 100000000 
                        + (uint64_t)(num[1] - '0')  * 1000000000
                        + (uint64_t)(num[2] - '0')  * 10000000000
                        + (uint64_t)(100000000000 * qtr[(num[1] - '0')* 10 + (num[2] - '0')]);                        

        if ((hashentry2 > (uint64_t)3999999999999))
        {                           
                printf("Bad Hash2 Value: %llu, %s\r\n",hashentry2,num);
                exit(0);
        }

        uint64_t hashentry3 =  ((((uint64_t)(num[3] - '0')))) 
                        + (uint64_t)(num[5] - '0')  * 10 
                        + (uint64_t)(num[7] - '0')  * 100 
                        + (uint64_t)(num[9] - '0')  * 1000
                        + (uint64_t)(num[11] - '0')   * 10000
                        + (uint64_t)(num[15] - '0')   * 100000
                        + (uint64_t)(num[17] - '0')   * 1000000
                        + (uint64_t)(num[9] - '0')   * 10000000                // 100M
                        + (uint64_t)(num[1] - '0')  * 100000000 
                        + (uint64_t)(num[4] - '0')  * 1000000000
                        + (uint64_t)(num[6] - '0')  * 10000000000
                        + (uint64_t)(100000000000 * qtr[(num[8] - '0')* 10 + (num[10] - '0')]);
            
                        

        if ((hashentry3 > (uint64_t)3999999999999))
        {                           
                printf("Bad Hash3 Value: %llu, %s\r\n",hashentry3,num);
                exit(0);
        }


        uint64_t hashentry4 =  ((((uint64_t)(num[14] - '0')))) 
                        + (uint64_t)(num[16] - '0')  * 10 
                        + (uint64_t)(num[18] - '0')  * 100 
                        + (uint64_t)(num[2] - '0')  * 1000
                        + (uint64_t)(num[4] - '0')  * 10000
                        + (uint64_t)(num[1] - '0')  * 100000
                        + (uint64_t)(num[2] - '0') * 1000000
                        + (uint64_t)(num[8] - '0') * 10000000                // 100M
                        + (uint64_t)(num[10] - '0') * 100000000 
                        + (uint64_t)(num[12] - '0') * 1000000000
                        + (uint64_t)(num[13] - '0') * 10000000000
                        + (uint64_t)(100000000000 * qtr[(num[14] - '0')* 10 + (num[15] - '0')]);
                
                        

        if ((hashentry4 > (uint64_t)3999999999999))
        {                           
                printf("Bad Hash4 Value: %llu, %s\r\n",hashentry4,num);
                exit(0);
        }

                BloomFlag1 = (S_TopBloom[hashentry1>>3] & (1<<((char)hashentry1&0x7)));
                BloomFlag2 = (S_TopBloom[hashentry2>>3] & (1<<((char)hashentry2&0x7)));
                BloomFlag3 = (S_TopBloom[hashentry3>>3] & (1<<((char)hashentry3&0x7)));
                BloomFlag4 = (S_TopBloom[hashentry4>>3] & (1<<((char)hashentry4&0x7)));

        if ((BloomFlag1==0) || (BloomFlag2==0) || (BloomFlag3==0) || (BloomFlag4==0))
        {
                //S_TopBloom[hashentry] = 1;
                S_TopBloom[hashentry1>>3] |= (1<<((char)hashentry1&0x7));
                S_TopBloom[hashentry2>>3] |= (1<<((char)hashentry2&0x7));
                S_TopBloom[hashentry3>>3] |= (1<<((char)hashentry3&0x7));
                S_TopBloom[hashentry4>>3] |= (1<<((char)hashentry4&0x7));
                S_UnallocatedHash = S_UnallocatedHash - 4;
        }
        else
        {       
                fprintf(datafile,"%s:%lld\r\n",num,offset);  
                S_BloomPass++;        
        }
}


void Complete()
{
        printf("Search Complete\n\r");
        gettimeofday(&S_stop, NULL);
        long timems = ((S_stop.tv_sec - S_start.tv_sec) * 1000) + ((S_stop.tv_usec - S_start.tv_usec)/1000);
        printf("Sequence Length %ld\r\n",S_digitsize);
        printf("Current Element Count: %6.2f M   - ",S_TotalElements/1000000.0);
        printf("Total Elements written to pass: %lld (%2.5f)\r\n",S_BloomPass,((float)S_BloomPass/(float)S_TotalElements));
        printf("Execution took %lu ms\n", timems);
        Cleanup();
        exit(0);
}

long totaltimer()
{
        struct timeval current;
        gettimeofday(&current, NULL);
        return (((current.tv_sec - S_start.tv_sec) * 1000) + ((current.tv_usec - S_start.tv_usec)/1000));
}

long lasttimer()
{
        struct timeval current;
        gettimeofday(&current, NULL);
        long r = ((current.tv_sec - S_last.tv_sec) * 1000) + ((current.tv_usec - S_last.tv_usec)/1000);
        S_last = current;
        return r;
}


