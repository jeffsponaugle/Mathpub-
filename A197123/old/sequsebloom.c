#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

// 
// A19 (1350168131352524443) occurs at 8858170606 and 4750253204 digits, 
// A20 (84756845106452435773) occurs at 17,601,613,331 and 15,494,062,638 digits.
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
#define HASHSIZE 1000000000             // Size of hashtable (1 billion)
#define PISOURCE "../pi-10billion.txt"
#define READBUFFERSIZE 41943040

// A197123
// The core structure used for each candidate number.
typedef struct s_candidate 
{
        char num[CANSIZE];              // The number itself (in ASCII)
        uint64_t offset;               // The offset into pi (digit location)
        struct s_candidate *next;       // Pointer to another candidate.
} candidate;

// Global variables ( used for performance )
candidate* *S_TopPreHash;               // Pre Hash table that holds bloom matches
uint64_t S_TotalElements = 0;
uint64_t S_TotalElementsLast = 0;          //  Total number of elements added to array.
uint64_t S_UnallocatedHash = 0;
long S_digitsize;                       //  The size 'digit match' we are looking for
struct timeval S_stop, S_start, S_last;         //  Start and stop time metrics.
long S_MaxDepth=0;
long S_MaxDepthCount=0;
char *S_NewBlockPtr=0;
long S_RemainingElementsinBlock=0;
char *S_NextBlockPtr;

FILE *logfile;
FILE *pisrc;
FILE *bloomsrc;
char S_digitbuffer[READBUFFERSIZE];            //(1024 x 4096)
int S_digitbuffercount = 0;
char *bloomfilename;

void Complete();
candidate* CreateCandidate(char* srcnum,uint64_t offset);
int CompareCandidates(candidate *a, candidate *b);
int CompareCandidatesRawAndUpdate(char* num, uint64_t numoffset, candidate *b);
void AddCandidate(candidate* *TopofHash, candidate* ca);
int TestAndUpdateCandidate(candidate* *TopofHash,char* num, uint64_t offset);
void PrintCandidateListSize();
long totaltimer();
long lasttimer();
candidate* GetElementMemory();
char getdigitbuffer();
void printnum(char* num);


int main(int argc, char *argv[])
{
        
        long digitfilter = 0;           // If we are filtering, this is the starting digit.
        int filterenabled = 0;          // 1 == filtering enabled
        char linebuffer[256];
        int bloomsrclines = 0;
        uint64_t offset = 0;
        char num[CANSIZE];

        logfile = fopen("pilog.txt","a");
        printf("\r\n\r\nStarting PI Search using bloom file.\r\n");
        char *p;
        if (argc==5)
        {
           digitfilter = strtol(argv[4], &p, 10); 
           printf("Limiting to [%ld]\r\n",digitfilter);
           filterenabled=1;
        }

        // Read input parameters
        // sequsebloom <bloomfile> <DIGITSTOTESTinM> <SIZEDIGITTOLOOKFOR> {FILTERNUMBER}
        
        // counter is the number of digits of PI we are going to search.
        uint64_t counter = (uint64_t)strtol(argv[2], &p, 10) * (1000000);

        // S_digitsize is the number of digits we are searching for.

        S_digitsize = strtol(argv[3], &p, 10);
        printf("Size to search for:%ld\r\n",S_digitsize);
        for (int u=0; u<CANSIZE;u++)
            num[u]=0; 

        printf("Processing %3.2f M numbers. \r\n",(float)counter/1000000.0);
        printf("Allocating pre-hash table\r\n");
        // Allocate the hashtable.
        S_TopPreHash = malloc(((uint64_t)sizeof(candidate*))*(uint64_t)HASHSIZE);
        printf("Allocated pre-hash table of size: %llu\r\n",(uint64_t)sizeof(candidate*)*HASHSIZE );
        printf("Buffer location: %llx\r\n",(uint64_t)S_TopPreHash);
        printf("Clearing hashtable space.\r\n");
        // Clear out all hash table entries
        memset(S_TopPreHash,0,HASHSIZE);
        S_UnallocatedHash = HASHSIZE;
        gettimeofday(&S_start, NULL);
        bloomfilename = argv[1];
        printf("Opening Bloom file:%s\r\n",bloomfilename);
        printf("Inserting pairs from Bloom File into Hash Array\r\n");
        bloomsrc = fopen(bloomfilename,"r");
        memset(linebuffer,0,128);
        while (fgets(linebuffer,128,bloomsrc) != NULL)
        {
            bloomsrclines++;
            if ( (strlen(linebuffer)>S_digitsize+4) & ( linebuffer[S_digitsize] == ':'))
            {    
                for (int i=0;i<S_digitsize;i++)
                     num[i] = (linebuffer[i] - '0');
                candidate* c = CreateCandidate(num,0);
                AddCandidate(S_TopPreHash,c);
            }
            else
            {
                printf("Error Parsing line: %s\r\n",linebuffer);
                exit(1);
            }
        }
        printf("\r\nTotal Bloom Src Lines: %d\r\n",bloomsrclines);
        //printf("last line: %s",linebuffer);
        //printf("parsed %s:%llu\r\n",num,offset);
        float hashratio = ((float)S_UnallocatedHash / (float) HASHSIZE)*100.0;
        printf("\r\nCurrent Element Count: %6.2fM : %5.2fM (Depth %ld Count %ld Unallocated %lld(%2.4f%%))   - \r\n",S_TotalElements/1000000.0,(float)(bloomsrclines/1000000.0),S_MaxDepth,S_MaxDepthCount,S_UnallocatedHash,hashratio);
        
        // clear the testing number.
        for (int u=0; u<CANSIZE;u++)
                num[u]=0; 
        S_TotalElements = 0;
        S_TotalElementsLast = 0;      
        S_MaxDepth=0;
        S_MaxDepthCount=0;

        printf("Starting search in Pi for match in Bloom Hash Table.\r\n");
    
        // ppi-10billion.txt contains the first 10 billion digits of pi, starting with
        //  3.14159265358.... 
        // We will skip over the '3.'

        pisrc = fopen(PISOURCE,"rb");
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
                        if (num[0] == ((char)digitfilter))
                        {
			   TestAndUpdateCandidate(S_TopPreHash,num,tries+1);
                        }
                }
                else
                {
       
	            TestAndUpdateCandidate(S_TopPreHash,num,tries+1);
                }

                // At this point the number has beed added, so if it makes it here it was not matched.
                
                // Every 1% of total tries we print a status report

                if (tries%(counter/100)==0)
                {
                         float hashratio = ((float)S_UnallocatedHash / (float) HASHSIZE)*100.0;
                         printf("\r\nCurrent Element Count: %6.2fK in Bloom out of : %5.2fK Digits (Depth %ld Count %ld Unallocated %llu(%2.4f%%))   - \r\n",S_TotalElements/1000.0,(float)(tries/1000.0),S_MaxDepth,S_MaxDepthCount,S_UnallocatedHash,hashratio);
                         long diff = totaltimer()+1;  // diff in ms
                         long difflast = lasttimer()+1;
                         float rate = S_TotalElements / (uint64_t)diff;  // rate in elem/ms, or kelem/s
                         float ratelast = (S_TotalElements-S_TotalElementsLast) / (uint64_t)difflast;
                         printf("       Cumulative Rate: %4.4f K Elements matched/sec   Section Rate: %4.4f K Elements matched/sec\r\n",rate,ratelast);
                         S_TotalElementsLast = S_TotalElements;

                        //getrusage(RUSAGE_SELF,&r_usage);
                        // Print the maximum resident set size used (in kilobytes).
                        //printf("Memory usage: %6.2f MB\n",r_usage.ru_maxrss/1000000.0);
                }

                // get the next digit, then move all the digits over by 1.

                nextdigit = getdigitbuffer();
                for (int t=1; t<S_digitsize; t++)
                        num[t-1]=num[t];
                num[S_digitsize-1] = nextdigit;
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
                fread(S_digitbuffer,READBUFFERSIZE,1,pisrc);
                S_digitbuffercount = READBUFFERSIZE;
        }
        return ((S_digitbuffer[READBUFFERSIZE-S_digitbuffercount--]) - '0');
}
candidate* GetElementMemory()
{
        char *r;

        if (S_RemainingElementsinBlock==0)
        {
                printf("*");
                S_NewBlockPtr = malloc(sizeof(candidate)*1024*1024);
                S_RemainingElementsinBlock = 1024*1024;
                S_NextBlockPtr = S_NewBlockPtr;
        }
        r = S_NextBlockPtr;
        S_NextBlockPtr = S_NextBlockPtr + sizeof(candidate);
        S_RemainingElementsinBlock--;
        return (candidate*)r;
}

candidate* CreateCandidate(char* srcnum,uint64_t offset)
{
        candidate *ca = GetElementMemory();
        ca->next=NULL;  
        ca->offset = offset;
        memcpy(&(ca->num),srcnum,CANSIZE);
        S_TotalElements++;
        return ca;
}

int CompareCandidates(candidate *a, candidate *b)
{
        char elementa, elementb;
        int offset = 0;

        elementa = a->num[offset];
        elementb = b->num[offset];
        while (elementa == elementb)
        {
                offset++;
                if (offset != S_digitsize)
                {
                        elementa = a->num[offset];
                        elementb = b->num[offset];
                }
                else
                {       
                        break;
                }

        }

        if (elementa > elementb) return 1;
        if (elementa < elementb) return -1;
        return 0;
}


int CompareCandidatesRawAndUpdate(char* num, uint64_t numoffset, candidate *b)
{
        char elementa, elementb;
        int index = 0;

        elementa = num[index];
        elementb = b->num[index];
        while (elementa == elementb)
        {
                index++;
                if (index != S_digitsize)
                {
                        elementa = num[index];
                        elementb = b->num[index];
                }
                else
                {       
                        break;
                }

        }

        if (elementa > elementb) return 1;
        if (elementa < elementb) return -1;
        // At this point we have found a matching number.  We will check to see if the 
        // the stored offset is 0.   If it is we will update the offset.
        if (b->offset == 0)
        {
                b->offset = numoffset;
		S_TotalElements++;
        }
        else 
        {
                // Winner!  If there is a non zero offset we much have a match.
                printf("Match found! \r\n");
                printnum(num);
                printf("\r\n");
                printnum(b->num);
                printf("  Locations: %lld,%lld\r\n",numoffset, b->offset);
                Complete();

        }
        return 0;
}

void printnum(char* num)
{
        for (int g=0;g<S_digitsize;g++)
        {
                printf("%c",num[g]+'0');
        }
}


void AddCandidate(candidate* *TopofHash,candidate* ca)
{
        candidate *nextc,*prevc, *S_Top;
        int r;
        long LocalDepth=1;

        // Lets start by figuring out which hash entry to use.    We used each digit in the 
        // number to pick a hash value.  The last value is used to select between the first 
        // and second 10B hash values. 

        uint64_t hashentry = (uint64_t)(ca->num[1]) 
                           + (uint64_t)(ca->num[2]) * 10 
                           + (uint64_t)(ca->num[3]) * 100 
                           + (uint64_t)(ca->num[4]) * 1000
                           + (uint64_t)(ca->num[5]) * 10000
                           + (uint64_t)(ca->num[6]) * 100000
                           + (uint64_t)(ca->num[7]) * 1000000
                           + (uint64_t)(ca->num[8]) * 10000000                // 10M
                           + (uint64_t)(ca->num[9]) * 100000000;
                          //+ (uint64_t)(ca->num[10] - '0') * 1000000000;


        if ((hashentry > (uint64_t)999999999))
        {
                printf("Bad Hash Value: %lld, %s\r\n",hashentry,ca->num);
                exit(0);
        }
                S_Top = TopofHash[hashentry];

        if (S_Top==0)
        {
                S_Top = ca;
                S_UnallocatedHash--;
        }
        else
        {       
                
                nextc=S_Top;
                prevc=S_Top;
                do {
                        r=CompareCandidates(ca,nextc);
                        if (r==0) 
                        {
                                printf("FIX: Match found! %s = %s\r\n",ca->num, nextc->num);
                                printf("Locations: %lld,%lld\r\n",ca->offset, nextc->offset);
                                Complete();
                        } else if (r==1)
                        {
                           // new element is greater than the current element, so lets skip to the next element.
                           prevc = nextc;
                           nextc = nextc->next; 
                           LocalDepth++;  

                        } else if (r==-1)
                        {
                           // new element is less than the current element, so we will add the element before the current one.
                           if (nextc==S_Top)
                           {
                                S_Top = ca;
                                ca->next = nextc;
                           } else
                           {
                                prevc->next = ca;
                                ca->next = nextc;
                           }  
                           LocalDepth++;
                           break;
                        }
                } while (nextc != 0);   // We will stop this loop once we reach a zero element
                // If nextc->next is zero, we must have gone through the entire loop, so lets just add to the end.
                if (nextc == 0)
                {
                        prevc->next = ca;
                        ca->next = 0;
                }
                if (S_MaxDepth < LocalDepth) 
                {
                        S_MaxDepth = LocalDepth;
                        S_MaxDepthCount = 1;
                }
                else if (S_MaxDepth == LocalDepth)
                {
                        S_MaxDepthCount++;
                }
                         
        }

        // Save the new top of this hash entry.
        TopofHash[hashentry] = S_Top;
}

int TestAndUpdateCandidate(candidate* *TopofHash,char* num, uint64_t offset)
{
        candidate *nextc,*prevc, *S_Top;
        int r;

        // Lets start by figuring out which hash entry to use.    We used each digit in the 
        // number to pick a hash value.  The last value is used to select between the first 
        // and second 10B hash values. 

        uint64_t hashentry = (uint64_t)(num[1]) 
                        + (uint64_t)(num[2]) * 10 
                        + (uint64_t)(num[3]) * 100 
                        + (uint64_t)(num[4]) * 1000
                        + (uint64_t)(num[5]) * 10000
                        + (uint64_t)(num[6]) * 100000
                        + (uint64_t)(num[7]) * 1000000
                        + (uint64_t)(num[8]) * 10000000                // 100M
                        + (uint64_t)(num[9]) * 100000000; 
                        //+ (uint64_t)(num[10] - '0') * 1000000000;

        if ((hashentry > (uint64_t)999999999) | (hashentry < (uint64_t)0))
        {
                printf("Bad Hash Value: %lld, %s\r\n",hashentry,num);
                exit(0);
        }

                S_Top = TopofHash[hashentry];

        if (S_Top==0)
        {   
                // If there is no entry for this hash value, there is not a match, so we
                // will ignore this add.
                return(0);
        }
        else
        {       
                
                nextc=S_Top;
                prevc=S_Top;
                do {
                        r=CompareCandidatesRawAndUpdate(num,offset,nextc);
                        if (r==0) 
                        {
                                return(1);
                        } else if (r==1)
                        {
                           // new element is greater than the current element, so lets skip to the next element.
                           prevc = nextc;
                           nextc = nextc->next;  

                        } else if (r==-1)
                        {
                           // new element is less than the current element so there is no match, return
                           return(0);
                        }
                } while (nextc != 0);   // We will stop this loop once we reach a zero element
                // If nextc->next is zero, we must have gone through the entire loop, so lets just add to the end.
                          
        }
        return(0);
}


void Complete()
{
        printf("Search Complete\n\r");
        gettimeofday(&S_stop, NULL);
        long timems = ((S_stop.tv_sec - S_start.tv_sec) * 1000) + ((S_stop.tv_usec - S_start.tv_usec)/1000);
        printf("Current Element Count: %6.2f M   - ",S_TotalElements/1000000.0);
        printf("Execution took %lu ms\n", timems);
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


