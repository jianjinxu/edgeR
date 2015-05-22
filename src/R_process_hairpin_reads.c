#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define BLOCKSIZE 10000000

typedef struct {
  char   *sequence;
  char   *sequence2;
  char   *sequenceRev;
  int    original_pos;
} a_barcode;

typedef struct {
   char   *sequence;
   int    original_pos;
} a_hairpin;

a_barcode **barcodes;
a_hairpin **hairpins;

int is_PairedReads; 
int is_DualIndexingReads;
int num_barcode;
int num_hairpin;
long num_read;
long **summary;
int barcode_start;
int barcode_end;
int barcode2_start;
int barcode2_end;
int barcode_start_rev;
int barcode_end_rev;
int barcode_length;
int barcode2_length;
int barcode_length_rev;
int hairpin_start;
int hairpin_end;
int hairpin_length;
int allow_shifting;
int shifting_n_base; 
int allow_mismatch;
int barcode_n_mismatch;
int hairpin_n_mismatch;
int allow_shifted_mismatch;
int isverbose;

long barcodecount;
long hairpincount;
long bchpcount;

int Get_Lines_In_File(FILE* fin) {
  int N=0, ch, last_ch='\n';
  while (ch=fgetc(fin)) { 
    if (ch=='\n') { ++N; }
    else if (ch==EOF) { 
      if (last_ch!='\n') { ++N; } // Capture non-newline-terminated last line.
      break;
    }
    last_ch=ch;
  }
  rewind(fin);
  return N;
}

void
Read_In_Barcodes(char* filename){
  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename,"r");
  
  // Getting number of lines in the file.
  num_barcode = Get_Lines_In_File(fin);
  barcodes=(a_barcode**)R_alloc(num_barcode+1, sizeof(a_barcode*));

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
    if (is_PairedReads > 0) {
      token = strtok(line, "\t");
      token = strtok(NULL, "\t");
      new_barcode->sequenceRev = (char *)malloc(barcode_length_rev * sizeof(char));
      strncpy(new_barcode->sequenceRev, token, barcode_length_rev);
    } else if (is_DualIndexingReads > 0) {
      token = strtok(line, "\t");
      token = strtok(NULL, "\t");
      new_barcode->sequence2 = (char *)malloc(barcode_length_rev * sizeof(char));
      strncpy(new_barcode->sequence2, token, barcode2_length);
    } else {
      new_barcode->sequenceRev = NULL;
    };
    barcodes[count] = new_barcode;
  }
  fclose(fin);
  free(line);

  Rprintf(" -- Number of Barcodes : %d\n", num_barcode);
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
 
  if (allow_mismatch > 0) {
    int i;
    for (i = 1; i <= num_barcode; i++){
      if (Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch) > 0) {
        return barcodes[i]->original_pos; 
      }   
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
  
  if (allow_mismatch > 0){
    int i;
    for (i = 1; i <= num_barcode; i++){
      if ((Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch) > 0) && 
          (Valid_Match(a_barcode_rev, barcodes[i]->sequenceRev, barcode_length_rev, barcode_n_mismatch) > 0)) {
         return barcodes[i]->original_pos;
      }
    }
  }

  return -1;
}


int
locate_barcode_dualIndexing(char *a_barcode, char *a_barcode2){
  int imin, imax, imid;
  imin = 1;
  imax = num_barcode;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    // compare forward barcode sequence 1 
    if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) > 0) {
      imax = imid - 1;
    } else {
      // same forward sequence 1, compare forward barcode sequence 2
      if (strncmp(barcodes[imid]->sequence2, a_barcode2, barcode2_length) < 0) {
        imin = imid + 1;
      } else if (strncmp(barcodes[imid]->sequence2, a_barcode2, barcode2_length) > 0) {
        imax = imid - 1;
      } else {
        return barcodes[imid]->original_pos;     
      } 
    }    
  }
  
  if (allow_mismatch > 0) {
    int i;
    for (i = 1; i <= num_barcode; i++){
      if ((Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch) > 0) && 
	    (Valid_Match(a_barcode2, barcodes[i]->sequence2, barcode2_length, barcode_n_mismatch) > 0)) {
        return barcodes[i]->original_pos;
      }
    }
  }
  
  return -1;
}

int
locate_exactmatch_hairpin(char *a_hairpin){
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
locate_mismatch_hairpin(char *a_hairpin){
  int i;
  for (i = 1; i <= num_hairpin; i++){
    if (Valid_Match(a_hairpin, hairpins[i]->sequence, hairpin_length, hairpin_n_mismatch) > 0) {
      return hairpins[i]->original_pos;
    }
  }
  return -1;
}


int 
locate_hairpin(char *a_hairpin, char *read){
  int hairpin_index;
  
  // check if a perfect match exists
  hairpin_index = locate_exactmatch_hairpin(a_hairpin); 
  if (hairpin_index > 0) {
    return hairpin_index;
  }
  
  // if a perfect match doesn't exist, check if a mismatch exists
  if (allow_mismatch > 0) {   
    hairpin_index = locate_mismatch_hairpin(a_hairpin);
	if (hairpin_index > 0) {
      return hairpin_index;
    }
  }

  if (allow_shifting > 0) {
    // Check if given hairpin can be mapped to a shifted location. 
    char *shiftedharipinstr;
    shiftedharipinstr = (char *)malloc(hairpin_length * sizeof(char));

    int index;
    // check shifting leftwards
    for (index = 1; index <= shifting_n_base; index++){
      strncpy(shiftedharipinstr, read + hairpin_start - 1 - index, hairpin_length);
      
	  hairpin_index = locate_exactmatch_hairpin(shiftedharipinstr);
	  
	  // check if mismatch on a shifted position is allowed and try to match
	  if ((hairpin_index <= 0) && (allow_shifted_mismatch)) {
	    hairpin_index = locate_mismatch_hairpin(shiftedharipinstr);
	  }
	  
      if (hairpin_index > 0) {
        return hairpin_index;
      } 
    }
    
	// check shifting rightwards
    for (index = 1; index <= shifting_n_base; index++){
      strncpy(shiftedharipinstr, read + hairpin_start - 1 + index, hairpin_length);
 
      hairpin_index = locate_exactmatch_hairpin(shiftedharipinstr); 
      
      if ((hairpin_index <= 0) && (allow_shifted_mismatch)) {	  
        hairpin_index = locate_mismatch_hairpin(shiftedharipinstr);
      } 

      if (hairpin_index > 0) {
        return hairpin_index;
      }
    }
  }
  
  return -1;
}

int 
barcode_compare(a_barcode *barcode1, a_barcode *barcode2){
  int ans;
  ans = strncmp(barcode1->sequence, barcode2->sequence, barcode_length);
  if (ans == 0) {
    if (is_PairedReads > 0){  
      ans = strncmp(barcode1->sequenceRev, barcode2->sequenceRev, barcode_length_rev);
    } else if (is_DualIndexingReads > 0){
      ans = strncmp(barcode1->sequence2, barcode2->sequence2, barcode2_length);
    }
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

  // Getting number of lines in the file.
  num_hairpin = Get_Lines_In_File(fin);
  hairpins=(a_hairpin**)R_alloc(num_hairpin+1, sizeof(a_hairpin*));

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
  free(line);

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
  char *this_barcode_2 = NULL;
  char *this_barcode_rev = NULL;
  char *this_hairpin = NULL;

  line = (char *)malloc(sizeof(char) * (len+1));
  fin = fopen(filename,"r");
  if (is_PairedReads > 0) {
    finRev = fopen(filename2, "r");
    line2 = (char *)malloc(sizeof(char) * (len+1));
  }

  if (isverbose > 0){
    if (is_PairedReads > 0) {
      Rprintf("Processing reads in %s and %s.\n", filename, filename2);
    } else {
      Rprintf("Processing reads in %s.\n", filename);
    }
  }

  this_barcode_for = (char *)malloc(barcode_length * sizeof(char));
  if (is_DualIndexingReads > 0) {
    this_barcode_2 = (char *)malloc(barcode2_length * sizeof(char));
  }
  if (is_PairedReads > 0) {
    this_barcode_rev = (char *)malloc(barcode_length_rev * sizeof(char));
  }
  this_hairpin = (char *)malloc(hairpin_length * sizeof(char));
  long line_count = 0;

  int barcode_index;
  int hairpin_index;

  while ((readline = fgets(line, len, fin)) != NULL){
    if (is_PairedReads > 0) {
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
    if (is_PairedReads > 0){    
      strncpy(this_barcode_rev, line2 + barcode_start_rev - 1, barcode_length_rev);
      barcode_index = locate_barcode_paired(this_barcode_for, this_barcode_rev); 
    } else if (is_DualIndexingReads > 0) {
      strncpy(this_barcode_2, line + barcode2_start - 1, barcode2_length);
      barcode_index = locate_barcode_dualIndexing(this_barcode_for, this_barcode_2);
    } else { 
      barcode_index = locate_barcode(this_barcode_for);
    }

    strncpy(this_hairpin, line + hairpin_start - 1, hairpin_length);
    hairpin_index = locate_hairpin(this_hairpin, line); 

    if (barcode_index > 0){
      barcodecount++;
    }

    if (hairpin_index > 0){
      hairpincount++;
    }
       
    if ((barcode_index > 0) && (hairpin_index > 0)) {
      summary[hairpin_index][barcode_index]++;
      bchpcount++;
    }

  } // end while

  if (isverbose > 0) {
    if (is_PairedReads > 0) {
      Rprintf("Number of reads in file %s and %s: %ld\n", filename, filename2, num_read_thisfile);  
    } else {
      Rprintf("Number of reads in file %s : %ld\n", filename, num_read_thisfile);  
    }
  }
 
  fclose(fin); 
  free(line);
  free(this_barcode_for);
  free(this_hairpin); 

  if (is_PairedReads > 0){  
    fclose(finRev);
    free(line2);
    free(this_barcode_rev);
  }
}


void 
Initialise(int IsPaired, int IsDualIndexing, 
           int barcodestart, int barcodeend, int barcode2start, int barcode2end, int barcodestartrev, int barcodeendrev, 
           int hairpinstart, int hairpinend, 
           int allowshifting, int shiftingnbase,
           int allowMismatch, int barcodemismatch, int hairpinmismatch, 
           int allowShiftedMismatch, int verbose){
			   
  num_barcode = 0;
  num_hairpin = 0;

  is_PairedReads = IsPaired;
  is_DualIndexingReads = IsDualIndexing;
  barcode_start = barcodestart;
  barcode_end = barcodeend;
  barcode2_start = barcode2start;
  barcode2_end = barcode2end;
  barcode_start_rev = barcodestartrev;
  barcode_end_rev = barcodeendrev;
  hairpin_start = hairpinstart;
  hairpin_end = hairpinend;
  barcode_length = barcode_end - barcode_start + 1;
  barcode2_length = barcode2_end - barcode2_start + 1;
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
    if (is_PairedReads > 0){
      free(barcodes[index]->sequenceRev);
    } 
    if (is_DualIndexingReads > 0){
      free(barcodes[index]->sequence2);
    }
    free(barcodes[index]);
  }
  for (index = 1; index <= num_hairpin; index++){
    free(hairpins[index]->sequence);
    free(hairpins[index]);
  }
  
  for (index = 0; index <= num_hairpin; index++){
    free(summary[index]);
  }
  free(summary);
}

void
Allocate_Summary_Table(void){
  int i, j;
  
  summary = (long **)malloc((num_hairpin+1) * sizeof(long *));
  for (i=0; i <= num_hairpin; i++){
    summary[i] = (long *)malloc((num_barcode+1) * sizeof(long));
  }
 
  for (i = 0; i <= num_hairpin; i++){
    for (j = 0; j <= num_barcode; j++){ 
      summary[i][j] = 0;  
	}
  }
}
void 
processHairpinReads(int *isPairedReads, int *isDualIndexingReads, 
                    char **file, char **file2, int *filecount, 
                    char**barcodeseqs, char**hairpinseqs, 
                    int *barcodestart, int *barcodeend, int *barcode2start, int *barcode2end, int *barcodestartrev, int *barcodeendrev, 
                    int *hairpinstart, int *hairpinend, 
                    int *allowShifting, int *shiftingnbase, 
                    int *allowMismatch, int *barcodemismatch, int *hairpinmismatch,
                    int *allowShiftedMismatch,
                    char **output, int *verbose)
{  
  Initialise(*isPairedReads, *isDualIndexingReads, 
             *barcodestart, *barcodeend, *barcode2start, *barcode2end, *barcodestartrev, *barcodeendrev, 
             *hairpinstart, *hairpinend, 
             *allowShifting, *shiftingnbase,
             *allowMismatch, *barcodemismatch, *hairpinmismatch, 
             *allowShiftedMismatch, *verbose);

  Read_In_Barcodes(*barcodeseqs);
  Sort_Barcodes(); // () is necessary to invoke function without parameters

  Read_In_Hairpins(*hairpinseqs);

  Check_Hairpins();
  Sort_Hairpins();  
  
  Allocate_Summary_Table();
  
  int i_file;

  for (i_file = 0; i_file < *filecount; i_file++){
    Process_Hairpin_Reads(file[i_file], file2[i_file]);
  }
 
  Rprintf("\nThe input run parameters are: \n");
  Rprintf(" -- Barcode: start position %d\t end position %d\t length %d\n", barcode_start, barcode_end, barcode_length);  
  if (is_DualIndexingReads){
    Rprintf(" -- Second Barcode in forward read: start position %d\t end position %d\t length %d\n", barcode2_start, barcode2_end, barcode2_length); 
  }
  if (is_PairedReads){
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
