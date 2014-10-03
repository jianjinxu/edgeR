#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define MAX_BARCODE 1000
#define MAX_HAIRPIN 100000
#define BLOCKSIZE 10000000

typedef struct {
   char   *sequence;
   char   *sequenceRev;
   int    original_pos;
} a_barcode;

typedef struct {
   char   *sequence;
   int    original_pos;
} a_hairpin;

a_barcode *barcodes[MAX_BARCODE];
a_hairpin *hairpins[MAX_HAIRPIN];

int isPairedReads;
int num_barcode;
int num_hairpin;
long num_read;
long hairpinreadcount[MAX_HAIRPIN];
long summary[MAX_HAIRPIN][MAX_BARCODE];
int barcode_start;
int barcode_end;
int barcode_start_rev;
int barcode_end_rev;
int barcode_length;
int barcode_length_rev;
int hairpin_start;
int hairpin_end;
int hairpin_length;
int allow_shifting;
int shifting_n_base; 
int allow_mismatch;
int num_mismatch_hairpin;
int barcode_n_mismatch;
int hairpin_n_mismatch;
int allow_shifted_mismatch;
int isverbose;

long barcodecount;
long hairpincount;
long bchpcount;

a_hairpin *mismatch_hairpins[MAX_HAIRPIN];

void
Read_In_Barcodes(char* filename){
  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename,"r");
  line = (char *)malloc(len+1);
  a_barcode *new_barcode;

  int count = 0;
  char * token;

  while ((readline = fgets(line, len, fin)) != NULL){
    count++;
    new_barcode = (a_barcode *)malloc(sizeof(a_barcode));
    new_barcode->sequence = (char *)malloc(barcode_length * sizeof(char));
    strncpy(new_barcode->sequence, line, barcode_length);
    new_barcode->original_pos = count;
    if (isPairedReads > 0) {
      token = strtok(line, "\t");
      token = strtok(NULL, "\t");
      new_barcode->sequenceRev = (char *)malloc(barcode_length_rev * sizeof(char));
      strncpy(new_barcode->sequenceRev, token, barcode_length_rev);
    } else {
      new_barcode->sequenceRev = NULL;
    };
    barcodes[count] = new_barcode;
  }
  fclose(fin);
  num_barcode = count;
  free(line);

  Rprintf(" -- Number of Barcodes : %d\n", num_barcode);
}

int
locate_barcode(char *a_barcode){
  int imin, imax, imid;
  imin = 1;
  imax = num_barcode;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    
    if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) > 0) {
      imax = imid - 1;
    } else {
      return barcodes[imid]->original_pos;
    }
  }
  return -1;  
}

int
locate_barcode_paired(char *a_barcode, char *a_barcode_rev){
  int imin, imax, imid;
  imin = 1;
  imax = num_barcode;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    // compare forward barcode sequence
    if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) > 0) {
      imax = imid - 1;
    } else {
      // same forward sequence, compare reverse barcode sequence
      if (strncmp(barcodes[imid]->sequenceRev, a_barcode_rev, barcode_length_rev) < 0) {
        imin = imid + 1;
      } else if (strncmp(barcodes[imid]->sequenceRev, a_barcode_rev, barcode_length_rev) > 0) {
        imax = imid - 1;
      } else {
        return barcodes[imid]->original_pos;     
      } 
    }    
  }
  return -1; 
}


int
locate_hairpin_impl(char *a_hairpin){
  int imin, imax, imid;
  imin = 1;
  imax = num_hairpin;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    
    if (strncmp(hairpins[imid]->sequence, a_hairpin, hairpin_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(hairpins[imid]->sequence, a_hairpin, hairpin_length) > 0) {
      imax = imid - 1;
    } else {
      return hairpins[imid]->original_pos;
    }
  }
  return -1;  
}

int
Valid_Match(char *sequence1, char *sequence2, int length, int threshold){
  int i_base;
  int mismatchbasecount = 0;
  for (i_base = 0; i_base < length; i_base++) {
    if (sequence1[i_base] != sequence2[i_base]) {
      mismatchbasecount++;		  
    }
  }
  if (mismatchbasecount <= threshold) {
    return 1;
  } else {
    return -1;
  }
}

int
locate_mismatch_barcode(char *a_barcode){
  int i;
  int match_index = -1;
  for (i = 1; i <= num_barcode; i++){
    if (Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch) > 0) {
      match_index = barcodes[i]->original_pos;
      break;
    }
  }
  return match_index;
}

int
locate_mismatch_barcode_paired(char *a_barcode, char *a_barcode_rev){
  int i;
  int match_index = -1;
  for (i = 1; i <= num_barcode; i++){
    if ((Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch) > 0) && 
	(Valid_Match(a_barcode_rev, barcodes[i]->sequenceRev, barcode_length_rev, barcode_n_mismatch) > 0)) {
      match_index = barcodes[i]->original_pos;
      break;
    }
  }
  return match_index;
}


int
locate_mismatch_hairpin(char *a_hairpin){
  int i;
  int match_index = -1;
  for (i = 1; i <= num_mismatch_hairpin; i++){
    if (Valid_Match(a_hairpin, mismatch_hairpins[i]->sequence, hairpin_length, hairpin_n_mismatch) > 0) {
      match_index = mismatch_hairpins[i]->original_pos;
      break;
    }
  }
  return match_index;
}


int 
locate_hairpin(char *a_hairpin, char *read, int doMismatch){
  int hairpin_index;
  if (doMismatch > 0){
    hairpin_index = locate_mismatch_hairpin(a_hairpin);
  } else {
    hairpin_index = locate_hairpin_impl(a_hairpin);
  }

  if ((hairpin_index <= 0) && (allow_shifting > 0)) {
    if ((doMismatch > 0) && (allow_shifted_mismatch <= 0)) {
      return hairpin_index;
      exit(0); 
    }
    // Check if given hairpin can be mapped to a shifted location. 
    char *shiftedharipinstr;
    shiftedharipinstr = (char *)malloc(hairpin_length * sizeof(char));

    int index;
    // check shifting leftwards
    for (index = 1; index <= shifting_n_base; index++){
      strncpy(shiftedharipinstr, read + hairpin_start - 1 - index, hairpin_length);
      if (doMismatch > 0){
        hairpin_index = locate_mismatch_hairpin(shiftedharipinstr);
      } else {
        hairpin_index = locate_hairpin_impl(shiftedharipinstr);
      }
      if (hairpin_index > 0) {
        break;
      }
    }
    // check shifting rightwards
    if (hairpin_index <= 0){
      for (index = 1; index <= shifting_n_base; index++){
        strncpy(shiftedharipinstr, read + hairpin_start - 1 + index, hairpin_length);
        
        if (doMismatch > 0){
          hairpin_index = locate_mismatch_hairpin(shiftedharipinstr);
        } else {
          hairpin_index = locate_hairpin_impl(shiftedharipinstr);
        }

        if (hairpin_index > 0) {
          break;
        }
      }
    }
  } 
  return hairpin_index;
}

int 
barcode_compare(a_barcode *barcode1, a_barcode *barcode2){
  int ans;
  ans = strncmp(barcode1->sequence, barcode2->sequence, barcode_length);
  if ((ans == 0) && (isPairedReads > 0)){  
    ans = strncmp(barcode1->sequenceRev, barcode2->sequenceRev, barcode_length_rev);
  }

  return ans;
}

void
Sort_Barcodes(void){
  int i, j;
  a_barcode *temp;
  for(i = 1; i < num_barcode; i++){
    for(j = i+1; j <= num_barcode; j++){
      if (barcode_compare(barcodes[i], barcodes[j]) > 0) {
	temp = barcodes[i];
	barcodes[i] = barcodes[j];
	barcodes[j] = temp;
      }
    }
  }
}

void
Read_In_Hairpins(char *filename){
  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename,"r");
  line = (char *)malloc(len+1);
  a_hairpin *new_hairpin;

  int count = 0;
  while ((readline = fgets(line, len, fin)) != NULL){
    count++;
    new_hairpin = (a_hairpin *)malloc(sizeof(a_hairpin));
    new_hairpin->sequence = (char *)malloc(hairpin_length * sizeof(char));
    new_hairpin->original_pos = count;
    strncpy(new_hairpin->sequence, line, hairpin_length);
    hairpins[count] = new_hairpin;
  }
  fclose(fin);
  num_hairpin = count;
  free(line);

  int i_hairpin;
  for (i_hairpin = 1; i_hairpin <= num_hairpin; i_hairpin++){
    hairpinreadcount[i_hairpin] = 0;
  }
  Rprintf(" -- Number of Hairpins : %d\n", num_hairpin);
}

void
Sort_Hairpins(void){
  int i, j;
  a_hairpin *temp;
  for(i = 1; i < num_hairpin; i++){
    for(j = i+1; j <= num_hairpin; j++){
      if (strncmp(hairpins[i]->sequence, hairpins[j]->sequence, hairpin_length) > 0){
	temp = hairpins[i];
	hairpins[i] = hairpins[j];
	hairpins[j] = temp;
      }
    }
  }
}

void
Process_Hairpin_Reads(char *filename, char *filename2){
  FILE *fin = NULL;
  FILE *finRev = NULL; 
  char * line = NULL;
  char * line2 = NULL;
  size_t len = 1000;
  char *readline;
  char *readline2;
  long num_read_thisfile = 0;

  char *this_barcode_for = NULL;
  char *this_barcode_rev = NULL;
  char *this_hairpin = NULL;

  line = (char *)malloc(sizeof(char) * (len+1));
  fin = fopen(filename,"r");
  if (isPairedReads > 0) {
    finRev = fopen(filename2, "r");
    line2 = (char *)malloc(sizeof(char) * (len+1));
  }

  if (isverbose > 0){
    if (isPairedReads > 0) {
      Rprintf("Processing reads in %s and %s.\n", filename, filename2);
    } else {
      Rprintf("Processing reads in %s.\n", filename);
    }
  }

  this_barcode_for = (char *)malloc(barcode_length * sizeof(char));
  if (isPairedReads > 0) {
    this_barcode_rev = (char *)malloc(barcode_length_rev * sizeof(char));
  }
  this_hairpin = (char *)malloc(hairpin_length * sizeof(char));
  long line_count = 0;

  int barcode_index;
  int hairpin_index;

  while ((readline = fgets(line, len, fin)) != NULL){
    if (isPairedReads > 0) {
      readline2 = fgets(line2, len, finRev);
	  if (readline2 == NULL) {
	    break;
	  }
    }
    line_count++;  
    if ((line_count % 4) != 2) {
      continue;
    }

    if ((isverbose > 0) && (num_read_thisfile % BLOCKSIZE == 0)) {
      Rprintf(" -- Processing %d million reads\n", (num_read_thisfile / BLOCKSIZE + 1) * 10);
    }
    num_read++;
    num_read_thisfile++;
  
    strncpy(this_barcode_for, line + barcode_start - 1, barcode_length);
    if (isPairedReads > 0){    
      strncpy(this_barcode_rev, line2 + barcode_start_rev - 1, barcode_length_rev);
      barcode_index = locate_barcode_paired(this_barcode_for, this_barcode_rev); 
    } else { 
      barcode_index = locate_barcode(this_barcode_for);
    }

    strncpy(this_hairpin, line + hairpin_start - 1, hairpin_length);
    hairpin_index = locate_hairpin(this_hairpin, line, 0); // not allowing mismatch

    if (barcode_index > 0){
      barcodecount++;
    }

    if (hairpin_index > 0){
      hairpincount++;
      hairpinreadcount[hairpin_index]++;
    }
       
    if ((barcode_index > 0) && (hairpin_index > 0)) {
      summary[hairpin_index][barcode_index]++;
      bchpcount++;
    }

  } // end while

  if (isverbose > 0) {
    if (isPairedReads > 0) {
      Rprintf("Number of reads in file %s and %s: %ld\n", filename, filename2, num_read_thisfile);  
    } else {
      Rprintf("Number of reads in file %s : %ld\n", filename, num_read_thisfile);  
    }
  }
 
  fclose(fin); 
  free(line);
  free(this_barcode_for);
  free(this_hairpin); 

  if (isPairedReads > 0){  
    fclose(finRev);
    free(line2);
    free(this_barcode_rev);
  }
}

void
Create_Mismatch_Hairpins_List(void){
  int i;
  num_mismatch_hairpin = 0;

  for (i = 1; i <= num_hairpin; i++){
    if (hairpinreadcount[i] == 0){
      num_mismatch_hairpin++;
      mismatch_hairpins[num_mismatch_hairpin] = hairpins[i];
    }
  }
  Rprintf("\nThere are %d hairpins without exact sequence match.\n", num_mismatch_hairpin);
}


void
Process_Mismatch(char *filename, char *filename2){

  FILE *fin = NULL;
  FILE *finRev = NULL;
  char * line = NULL;
  char * line2 = NULL;
  size_t len = 1000;
  char *readline;
  char *readline2;
  long num_read_thisfile = 0;

  char * this_hairpin;
  char * this_barcode_for;
  char * this_barcode_rev = NULL;

  line = (char *)malloc(len+1);
  fin = fopen(filename,"r");
  if (isPairedReads > 0) {
    finRev = fopen(filename2, "r");
    line2 = (char *)malloc(sizeof(char) * (len+1));
  }
  if (isverbose > 0){
    if (isPairedReads > 0) {
      Rprintf("Re-processing reads in %s and %s, considering sequence mismatch\n", filename, filename2);
    } else {
      Rprintf("Re-processing reads in %s, considering sequence mismatch\n", filename);
    }
  }

  this_barcode_for = (char *)malloc(barcode_length * sizeof(char));
  if (isPairedReads > 0) {
    this_barcode_rev = (char *)malloc(barcode_length_rev * sizeof(char));
  }
  this_hairpin = (char *)malloc(hairpin_length * sizeof(char));

  long line_count = 0;

  int barcode_index;
  int hairpin_index;

  int new_barcode_index;
  int new_hairpin_index;

  while ((readline = fgets(line, len, fin)) != NULL){
    if (isPairedReads > 0) {
      readline2 = fgets(line2, len, finRev);
	  if (readline2 == NULL) {
	    break;
	  }
    }

    line_count++;
    if ((line_count % 4) != 2) {
      continue;
    }

    if ((isverbose > 0) && (num_read_thisfile % BLOCKSIZE == 0)) {
      Rprintf(" -- Processing %d million reads\n", (num_read_thisfile / BLOCKSIZE + 1) * 10);
    }
    num_read++;
    num_read_thisfile++;

    // re-do the mapping in Process_Hairpin_Reads()
    strncpy(this_barcode_for, line + barcode_start - 1, barcode_length);
    if (isPairedReads > 0){    
      strncpy(this_barcode_rev, line2 + barcode_start_rev - 1, barcode_length_rev);
      barcode_index = locate_barcode_paired(this_barcode_for, this_barcode_rev); 
    } else { 
      barcode_index = locate_barcode(this_barcode_for);
    }

    strncpy(this_hairpin, line + hairpin_start - 1, hairpin_length);
    hairpin_index = locate_hairpin(this_hairpin, line, 0); //not allowing mismatch

    // only re-process reads without perfect hairpin match or without perfect barcode match;    
    if ((barcode_index > 0) && (hairpin_index > 0)) {
      continue;
    }

    if (barcode_index > 0) {
      new_barcode_index = barcode_index;
    } else {
      if (isPairedReads > 0){
        new_barcode_index = locate_mismatch_barcode_paired(this_barcode_for, this_barcode_rev);
      } else {
        new_barcode_index = locate_mismatch_barcode(this_barcode_for);
      }

      if (new_barcode_index > 0) {
	barcodecount++;
      }
    } 
 
    // re-match hairpin:
    if (hairpin_index > 0){
      new_hairpin_index = hairpin_index;
    } else {
      new_hairpin_index = locate_hairpin(this_hairpin, line, 1);

      if (new_hairpin_index > 0) {
        hairpincount++;
      } 
    }
  
    if ((new_barcode_index > 0) && (new_hairpin_index > 0)) {
      summary[new_hairpin_index][new_barcode_index]++;
      bchpcount++;
    }
  } // end while 

  fclose(fin); 
  free(line);
  free(this_barcode_for);
  free(this_hairpin); 

  if (isPairedReads > 0){  
    fclose(finRev);
    free(line2);
    free(this_barcode_rev);
  }

}


void 
Initialise(int IsPaired, int barcodestart, int barcodeend, int barcodestartrev, int barcodeendrev, int hairpinstart, int hairpinend, 
	   int allowshifting, int shiftingnbase,
	   int allowMismatch, int barcodemismatch, int hairpinmismatch, 
	   int allowShiftedMismatch, int verbose){
  int i, j;
  for(i = 0; i < MAX_HAIRPIN; i++) {
    for(j = 0; j < MAX_BARCODE; j++) {
      summary[i][j] = 0;
    }
  }
  num_barcode = 0;
  num_hairpin = 0;

  isPairedReads = IsPaired;
  barcode_start = barcodestart;
  barcode_end = barcodeend;
  barcode_start_rev = barcodestartrev;
  barcode_end_rev = barcodeendrev;
  hairpin_start = hairpinstart;
  hairpin_end = hairpinend;
  barcode_length = barcode_end - barcode_start + 1;
  barcode_length_rev = barcode_end_rev - barcode_start_rev + 1;
  hairpin_length = hairpin_end - hairpin_start + 1;

  allow_shifting = allowshifting;
  shifting_n_base = shiftingnbase;
  allow_mismatch = allowMismatch;
  barcode_n_mismatch = barcodemismatch;
  hairpin_n_mismatch = hairpinmismatch;
  allow_shifted_mismatch = allowShiftedMismatch;
  isverbose = verbose;

  num_read = 0;
  barcodecount = 0;
  hairpincount = 0;
  bchpcount = 0;
}

void
Output_Summary_Table(char *output){
  int i, j;
  FILE *fout;
  fout = fopen(output, "w");
  for(i = 1; i <= num_hairpin; i++) {
    fprintf(fout, "%ld", summary[i][1]);
    for(j = 2; j <= num_barcode; j++) {
      fprintf(fout, "\t%ld", summary[i][j]);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}

void
Check_Hairpins(void){
  int p, q;
  char base;
  for(p = 1; p <= num_hairpin; p++){
    for(q = 0; q < hairpin_length; q++){
      base = hairpins[p]->sequence[q];
      if ((base != 'A') && (base != 'T') && (base != 'G') && (base != 'C')){
	Rprintf("Hairpin no.%d: %s contains invalid base %c\n", p, hairpins[p]->sequence, base);
      }
    }
  }
}

void
Clean_Up(void){
  int index;
  for (index = 1; index <= num_barcode; index++){
    free(barcodes[index]->sequence);
    if (isPairedReads > 0){
      free(barcodes[index]->sequenceRev);
    }
    free(barcodes[index]);
  }
  for (index = 1; index <= num_hairpin; index++){
    free(hairpins[index]->sequence);
    free(hairpins[index]);
  }
}


void 
processHairpinReads(int *isPairedReads,
                    char **file, char **file2, int *filecount, 
		    char**barcodeseqs, char**hairpinseqs, 
		    int *barcodestart, int *barcodeend, int *barcodestartrev, int *barcodeendrev, int *hairpinstart, int *hairpinend, 
		    int *allowShifting, int *shiftingnbase, 
		    int *allowMismatch, int *barcodemismatch, int *hairpinmismatch,
		    int *allowShiftedMismatch,
		    char **output, int *verbose)
{  
  Initialise(*isPairedReads, *barcodestart, *barcodeend, *barcodestartrev, *barcodeendrev, *hairpinstart, *hairpinend, 
	     *allowShifting, *shiftingnbase,
	     *allowMismatch, *barcodemismatch, *hairpinmismatch, 
	     *allowShiftedMismatch, *verbose);

  Read_In_Barcodes(*barcodeseqs);
  Sort_Barcodes(); // () is necessary to invoke function without parameters

  Read_In_Hairpins(*hairpinseqs);

  Check_Hairpins();
  Sort_Hairpins();  
  
  int i_file;

  for (i_file = 0; i_file < *filecount; i_file++){
    Process_Hairpin_Reads(file[i_file], file2[i_file]);
  }
 
  if (allow_mismatch > 0){
    num_read = 0; 
    // reset total number of read to 0, recheck initial barcode/hairpin index
    Create_Mismatch_Hairpins_List();
    if (num_mismatch_hairpin > 0){
      for (i_file = 0; i_file < *filecount; i_file++){
	if (isPairedReads > 0) {
          Process_Mismatch(file[i_file], file2[i_file]);
        } else {
          Process_Mismatch(file[i_file], NULL);
        }
      }
    }  
  }

  Rprintf("\nThe input run parameters are: \n");
  Rprintf(" -- Barcode: start position %d\t end position %d\t length %d\n", barcode_start, barcode_end, barcode_length);  
  if (isPairedReads > 0){
    Rprintf(" -- Barcode in reverse read: start position %d\t end position %d\t length %d\n", barcode_start_rev, barcode_end_rev, barcode_length_rev); 
  }
  Rprintf(" -- Hairpin: start position %d\t end position %d\t length %d\n", hairpin_start, hairpin_end, hairpin_length); 
  if (allow_shifting) {
    Rprintf(" -- Allow hairpin sequences to be matched to a shifted position, <= %d base left or right of the specified positions. \n", shifting_n_base);
  } else {
    Rprintf(" -- Hairpin sequences need to match at specified positions. \n");
  }
  if (allow_mismatch) {
    Rprintf(" -- Allow sequence mismatch, <= %d base in barcode sequence and <= %d base in hairpin sequence. \n", barcode_n_mismatch, hairpin_n_mismatch );
  } else {
    Rprintf(" -- Mismatch in barcode/hairpin sequences not allowed. \n");
  } 

  Rprintf("\nTotal number of read is %ld \n", num_read);
  Rprintf("There are %ld reads (%.4f percent) with barcode matches\n", barcodecount, 100.0*barcodecount/num_read);
  Rprintf("There are %ld reads (%.4f percent) with hairpin matches\n", hairpincount, 100.0*hairpincount/num_read);
  Rprintf("There are %ld reads (%.4f percent) with both barcode and hairpin matches\n", bchpcount, 100.0*bchpcount/num_read);

  Output_Summary_Table(*output);

  Clean_Up();
}






       
