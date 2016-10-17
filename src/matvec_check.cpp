#include "matvec_check.h"

matvec_check::matvec_check(SEXP incoming, int nr, int nc) : nrow(nr), ncol(nc), mycheck(NULL), index(0), libdex(0), temp(NULL), toget(NULL) { 

    SEXP dims=getAttrib(incoming, R_DimSymbol);
    if (!isInteger(dims) || LENGTH(dims)!=2) { 
        throw std::runtime_error("matrix dimensions should be an integer vector of length 2");
    }
    const int NR=INTEGER(dims)[0], NC=INTEGER(dims)[1];
    if (LENGTH(incoming)!=NC*NR) {
        throw std::runtime_error("recorded matrix dimensions are not consistent with its length"); 
    }

    SEXP rep_r=getAttrib(incoming, install("repeat.row"));
    if (!isLogical(rep_r) || LENGTH(rep_r)!=1) {
        throw std::runtime_error("repeat_row specification must be a logical scalar");
    }
    repeat_row=asLogical(rep_r);
    if (!repeat_row) {
        if (NR!=nrow){ 
            throw std::runtime_error("matrix dimensions are not consistent for non-repeating number of rows");
        }
    } else if (NR!=1) {
        throw std::runtime_error("only one row should be present if it is repeating");
    }
    maxdex=NR;

    SEXP rep_c=getAttrib(incoming, install("repeat.col"));
    if (!isLogical(rep_c) || LENGTH(rep_c)!=1) {
        throw std::runtime_error("repeat_col specification must be a logical scalar");
    }
    repeat_col=asLogical(rep_c);
    if (!repeat_col) {
        if (NC!=ncol){ 
            throw std::runtime_error("matrix dimensions are not consistent for non-repeating number of columns");
        }
    } else if (NC!=1) {
        throw std::runtime_error("only one column should be present if it is repeating");
    }

	if (!isReal(incoming)) {
		throw std::runtime_error("matrix should be double-precision"); 
	}	
	mycheck=REAL(incoming);
    temp=new double[ncol];
    
    try {
        // Setting the initial value of 'temp'.
        if (repeat_row) {
            if (repeat_col) {
                std::fill(temp, temp+ncol, *mycheck);
            } else {
                temp=new double[ncol];
                std::copy(mycheck, mycheck+ncol, temp);
            }
        } else {
            advance();
        }
    } catch (std::exception& e) {
        delete [] temp;
        throw;
    }
    
	return;
}

void matvec_check::advance() {
    if (repeat_row || index>=maxdex) {
        // No need for updating, if repeated across rows; or if we're at the end.
        return;
    }
    if (repeat_col) {  
        // Updating for new row.
        std::fill(temp, temp+ncol, mycheck[index]);
    } else {
        // Updating for new row/col.
        toget=mycheck + index;
        for (libdex=0; libdex<ncol; ++libdex) {
            temp[libdex]=*toget;
            toget+=nrow;
        }
    }
    ++index;
    return;
}

const double* const matvec_check::access() const { 
    return temp;
}

const bool matvec_check::is_row_repeated () const {
    return repeat_row;
}

const bool matvec_check::is_col_repeated () const {
    return repeat_col;
}

matvec_check::~matvec_check() {
	if (temp!=NULL) { delete [] temp; }
}

count_holder::count_holder (SEXP y) : yiptr(NULL), ydptr(NULL) {
    SEXP dims=getAttrib(y, R_DimSymbol);
    if (!isInteger(dims) || LENGTH(dims)!=2) {                         
        throw std::runtime_error("matrix dimensions should be an integer vector of length 2");                            
    }                
    num_tags=INTEGER(dims)[0];                    
    num_libs=INTEGER(dims)[1];

    is_integer=isInteger(y);                                        
    if (is_integer) {                                                             
        yiptr=INTEGER(y);                                                                     
    } else if (isReal(y)) {                                                                                        
        ydptr=REAL(y);
    } else {                                                                                                                    
        throw std::runtime_error("count matrix must be integer or double-precision");
    }

    if (LENGTH(y)!=num_tags*num_libs) {
        throw std::runtime_error("dimensions of the count matrix are not as specified");
    } 
    return;
}

void count_holder::fill_and_next(double* yptr) {
    if (is_integer) {
        for (libdex=0; libdex<num_libs; ++libdex) { yptr[libdex]=double(yiptr[libdex*num_tags]); }
        ++yiptr;
    } else {
        for (libdex=0; libdex<num_libs; ++libdex) { yptr[libdex]=ydptr[libdex*num_tags]; }
        ++ydptr;
    }
    return;
}

bool count_holder::is_data_integer() const { return is_integer; }
const int* count_holder::get_raw_int () const { return yiptr; }
const double* count_holder::get_raw_double () const { return ydptr; }

int count_holder::get_ntags() const { return num_tags; }
int count_holder::get_nlibs() const { return num_libs; }
