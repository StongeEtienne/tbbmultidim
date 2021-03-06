/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlFastOps.h"

#include <vcl_cstdlib.h> // abort()
#include <vcl_cstring.h> // memset()
#include <vcl_iostream.h>

namespace crl {
/*
void FastOps::ggt( vnl_matrix<double> &out, const vnl_vector<double>& g )
{
	double const* bb = g.data_block();
	int n = g.size();

	out.set_size(n,n);
	out.fill(0);
	for ( int i = 0; i < n; ++i)
		for ( int j = 0; j < n; ++j) {
			out[i][j] += bb[j] * bb[i];
		}	

}

void FastOps::AtBA(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
{
  const unsigned int na = A.columns();
  const unsigned int mb = B.rows();
  const unsigned int ma = A.rows();
  const unsigned int nb = B.columns();

  // Verify matrices compatible
  if (ma != mb) {
    vcl_cerr << "vnl_fastops::ABAt: argument sizes do not match: " << ma << " != " << mb << '\n';
    vcl_abort();
  }

  // Verify matrices compatible
  if (ma != nb) {
    vcl_cerr << "vnl_fastops::ABAt: argument sizes do not match: " << ma << " != " << nb << '\n';
    vcl_abort();
  }

  // Verify output is the right size
  if (out.rows() != na || out.columns() != na)
    out.set_size(na,na);

  double const* const* a = A.data_array();
  double const* const* b = B.data_array();
  double** outdata = out.data_array();

  // initialize
  for (unsigned int i = 0; i < ma; ++i)
    for (unsigned int w = 0; w < ma; ++w)
      outdata[i][w] = 0.0;

  for (unsigned int i = 0; i < na; ++i)
    for (unsigned int j = 0; j < nb; ++j) {
      double accum = 0;

      for (unsigned int k = 0; k < ma; ++k)
        accum += a[k][i] * b[k][j];
      for (unsigned int w = 0; w < na; ++w)
        outdata[i][w] += accum * a[j][w];
    }
}

void FastOps::AtBA_3x3(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
{

  // Verify output is the right size (there is a test inside set_size, so if 3x3, don t modify)
  out.set_size(3,3);

  double const* const* a = A.data_array();
  double const* const* b = B.data_array();
  double** outdata = out.data_array();

  // initialize
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int w = 0; w < 3; ++w)
      outdata[i][w] = 0.0;

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j) {
      double accum = 0;

      for (unsigned int k = 0; k < 3; ++k)
        accum += a[k][i] * b[k][j];
      for (unsigned int w = 0; w < 3; ++w)
        outdata[i][w] += accum * a[j][w];
    }
}



void FastOps::AtBA_3x3(vnl_matrix_ref<double> out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
{

  // Verify output is the right size (there is a test inside set_size, so if 3x3, don t modify)
  //out.set_size(3,3);

  double const* const* a = A.data_array();
  double const* const* b = B.data_array();
  double** outdata = out.data_array();

  // initialize
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int w = 0; w < 3; ++w)
      outdata[i][w] = 0.0;

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j) {
      double accum = 0;

      for (unsigned int k = 0; k < 3; ++k)
        accum += a[k][i] * b[k][j];
      for (unsigned int w = 0; w < 3; ++w)
        outdata[i][w] += accum * a[j][w];
    }
}

void FastOps::AtBA_Bsym_3x3(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
{
	// Verify output is the right size (there is a test inside set_size, so if 3x3, don t modify)
	out.set_size(3,3);

	double const* const* a = A.data_array();
	double const* const* b = B.data_array();
	double** outdata = out.data_array();

	outdata[0][0] = a[2][0]*(a[2][0]*b[2][2]+a[1][0]*b[1][2]+a[0][0]*b[0][2])+a[1][0]*(b[1][2]*a[2][0]+a[1][0]*b[1][1]+a[0][0]*b[0][1])+a[0][0]*(b[0][2]*a[2][0]+b[0][1]*a[1][0]+a[0][0]*b[0][0]);
	outdata[0][1] = a[2][0]*(a[2][1]*b[2][2]+a[1][1]*b[1][2]+a[0][1]*b[0][2])+a[1][0]*(b[1][2]*a[2][1]+a[1][1]*b[1][1]+a[0][1]*b[0][1])+a[0][0]*(b[0][2]*a[2][1]+b[0][1]*a[1][1]+b[0][0]*a[0][1]);
	outdata[0][2] = a[2][0]*(a[2][2]*b[2][2]+a[1][2]*b[1][2]+a[0][2]*b[0][2])+a[1][0]*(b[1][2]*a[2][2]+b[1][1]*a[1][2]+b[0][1]*a[0][2])+a[0][0]*(b[0][2]*a[2][2]+b[0][1]*a[1][2]+b[0][0]*a[0][2]);
	outdata[1][1] = a[2][1]*(a[2][1]*b[2][2]+a[1][1]*b[1][2]+a[0][1]*b[0][2])+a[1][1]*(b[1][2]*a[2][1]+a[1][1]*b[1][1]+a[0][1]*b[0][1])+a[0][1]*(b[0][2]*a[2][1]+b[0][1]*a[1][1]+b[0][0]*a[0][1]);
	outdata[1][2] = a[2][1]*(a[2][2]*b[2][2]+a[1][2]*b[1][2]+a[0][2]*b[0][2])+a[1][1]*(b[1][2]*a[2][2]+b[1][1]*a[1][2]+b[0][1]*a[0][2])+a[0][1]*(b[0][2]*a[2][2]+b[0][1]*a[1][2]+b[0][0]*a[0][2]);
	outdata[2][2] = a[2][2]*(a[2][2]*b[2][2]+a[1][2]*b[1][2]+a[0][2]*b[0][2])+a[1][2]*(b[1][2]*a[2][2]+b[1][1]*a[1][2]+b[0][1]*a[0][2])+a[0][2]*(b[0][2]*a[2][2]+b[0][1]*a[1][2]+b[0][0]*a[0][2]) ;

	outdata[1][0] = outdata[0][1] ;
	outdata[2][0] = outdata[0][2] ; 
	outdata[2][1] = outdata[1][2] ;
}

void FastOps::AtBA_Bsym_3x3(vnl_matrix_ref<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
{
	double const* const* a = A.data_array();
	double const* const* b = B.data_array();
	double** outdata = out.data_array();

	outdata[0][0] = a[2][0]*(a[2][0]*b[2][2]+a[1][0]*b[1][2]+a[0][0]*b[0][2])+a[1][0]*(b[1][2]*a[2][0]+a[1][0]*b[1][1]+a[0][0]*b[0][1])+a[0][0]*(b[0][2]*a[2][0]+b[0][1]*a[1][0]+a[0][0]*b[0][0]);
	outdata[0][1] = a[2][0]*(a[2][1]*b[2][2]+a[1][1]*b[1][2]+a[0][1]*b[0][2])+a[1][0]*(b[1][2]*a[2][1]+a[1][1]*b[1][1]+a[0][1]*b[0][1])+a[0][0]*(b[0][2]*a[2][1]+b[0][1]*a[1][1]+b[0][0]*a[0][1]);
	outdata[0][2] = a[2][0]*(a[2][2]*b[2][2]+a[1][2]*b[1][2]+a[0][2]*b[0][2])+a[1][0]*(b[1][2]*a[2][2]+b[1][1]*a[1][2]+b[0][1]*a[0][2])+a[0][0]*(b[0][2]*a[2][2]+b[0][1]*a[1][2]+b[0][0]*a[0][2]);
	outdata[1][1] = a[2][1]*(a[2][1]*b[2][2]+a[1][1]*b[1][2]+a[0][1]*b[0][2])+a[1][1]*(b[1][2]*a[2][1]+a[1][1]*b[1][1]+a[0][1]*b[0][1])+a[0][1]*(b[0][2]*a[2][1]+b[0][1]*a[1][1]+b[0][0]*a[0][1]);
	outdata[1][2] = a[2][1]*(a[2][2]*b[2][2]+a[1][2]*b[1][2]+a[0][2]*b[0][2])+a[1][1]*(b[1][2]*a[2][2]+b[1][1]*a[1][2]+b[0][1]*a[0][2])+a[0][1]*(b[0][2]*a[2][2]+b[0][1]*a[1][2]+b[0][0]*a[0][2]);
	outdata[2][2] = a[2][2]*(a[2][2]*b[2][2]+a[1][2]*b[1][2]+a[0][2]*b[0][2])+a[1][2]*(b[1][2]*a[2][2]+b[1][1]*a[1][2]+b[0][1]*a[0][2])+a[0][2]*(b[0][2]*a[2][2]+b[0][1]*a[1][2]+b[0][0]*a[0][2]) ;

	outdata[1][0] = outdata[0][1] ;
	outdata[2][0] = outdata[0][2] ; 
	outdata[2][1] = outdata[1][2] ;
}

void FastOps::AtBA_Bdiag_3x3(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_diag_matrix<double>& B)
{
	// Verify output is the right size (there is a test inside set_size, so if 3x3, don t modify)
	out.set_size(3,3);

	double const* const* a = A.data_array();
	double const* b = B.data_block();
	double** outdata = out.data_array();

	outdata[0][0] = a[2][0]*a[2][0]*b[2]+a[1][0]*a[1][0]*b[1]+a[0][0]*a[0][0]*b[0];
	outdata[0][1] = a[2][0]*a[2][1]*b[2]+a[1][0]*a[1][1]*b[1]+a[0][0]*b[0]*a[0][1];
	outdata[0][2] = a[2][0]*a[2][2]*b[2]+a[1][0]*b[1]*a[1][2]+a[0][0]*b[0]*a[0][2];
	outdata[1][1] = a[2][1]*a[2][1]*b[2]+a[1][1]*a[1][1]*b[1]+b[0]*a[0][1]*a[0][1];
	outdata[1][2] = a[2][1]*a[2][2]*b[2]+a[1][1]*b[1]*a[1][2]+b[0]*a[0][1]*a[0][2];
	outdata[2][2] = a[2][2]*a[2][2]*b[2]+b[1]*a[1][2]*a[1][2]+b[0]*a[0][2]*a[0][2];

	outdata[1][0] = outdata[0][1] ;
	outdata[2][0] = outdata[0][2] ; 
	outdata[2][1] = outdata[1][2] ;
}

void FastOps::AtBA_Bdiag_3x3(vnl_matrix_fixed<double,3,3>& out, const vnl_matrix_fixed<double,3,3>& A, const vnl_vector_fixed<double,3>& B)
{
	// Verify output is the right size (there is a test inside set_size, so if 3x3, don t modify)

    double const* a = A.data_block();
    double const* b = B.data_block();
    double* outdata = out.data_block();

    double temp10 = a[6]*a[7]*b[2]+a[3]*a[4]*b[1]+a[0]*b[0]*a[1];
    double temp20 = a[6]*a[8]*b[2]+a[3]*b[1]*a[5]+a[0]*b[0]*a[2];
    double temp21 = a[7]*a[8]*b[2]+a[4]*b[1]*a[5]+b[0]*a[1]*a[2];

    outdata[1] = temp10;
    outdata[3] = temp10;
    outdata[2] = temp20;
    outdata[6] = temp20;
    outdata[7] = temp21;
    outdata[5] = temp21;

    outdata[0] = a[6]*a[6]*b[2]+a[3]*a[3]*b[1]+a[0]*a[0]*b[0];
    outdata[4] = a[7]*a[7]*b[2]+a[4]*a[4]*b[1]+b[0]*a[1]*a[1];
    outdata[8] = a[8]*a[8]*b[2]+b[1]*a[5]*a[5]+b[0]*a[2]*a[2];

}

void FastOps::AtBA_Bdiag_3x3(vnl_matrix_ref<double> out, const vnl_matrix<double>& A, const vnl_diag_matrix<double>& B)
{
	double const* const* a = A.data_array();
	double const* b = B.data_block();
	double** outdata = out.data_array();

	outdata[0][0] = a[2][0]*a[2][0]*b[2]+a[1][0]*a[1][0]*b[1]+a[0][0]*a[0][0]*b[0];
	outdata[0][1] = a[2][0]*a[2][1]*b[2]+a[1][0]*a[1][1]*b[1]+a[0][0]*b[0]*a[0][1];
	outdata[0][2] = a[2][0]*a[2][2]*b[2]+a[1][0]*b[1]*a[1][2]+a[0][0]*b[0]*a[0][2];
	outdata[1][1] = a[2][1]*a[2][1]*b[2]+a[1][1]*a[1][1]*b[1]+b[0]*a[0][1]*a[0][1];
	outdata[1][2] = a[2][1]*a[2][2]*b[2]+a[1][1]*b[1]*a[1][2]+b[0]*a[0][1]*a[0][2];
	outdata[2][2] = a[2][2]*a[2][2]*b[2]+b[1]*a[1][2]*a[1][2]+b[0]*a[0][2]*a[0][2];

	outdata[1][0] = outdata[0][1] ;
	outdata[2][0] = outdata[0][2] ; 
	outdata[2][1] = outdata[1][2] ;
}


double FastOps::btAb_Asym_3x3 (const vnl_matrix<double>& A, const vnl_vector<double>& b)
{
	double const* const* a = A.data_array();
	double const* v = b.data_block();

	//v[2]*(v[2]*a[2,2]+v[1]*a[1,2]+v[0]*a[0,2])+v[1]*(a[1,2]*v[2]+v[1]*A[1,1]+v[0]*a[0,1])+v[0]*(a[0,2]*v[2]+a[0,1]*v[1]+v[0]*a[0,0])
 return (
	  v[2]*(v[2]*a[2][2]+v[1]*a[1][2]+v[0]*a[0][2])
	+ v[1]*(a[1][2]*v[2]+v[1]*a[1][1]+v[0]*a[0][1])
	+ v[0]*(a[0][2]*v[2]+a[0][1]*v[1]+v[0]*a[0][0]) );
}

double FastOps::btAb_Asym_3x3 (const vnl_matrix<double>& A, const vnl_vector_fixed<double,3>& b)
{
    double const* const* a = A.data_array();
	double const* v = b.data_block();

	//v[2]*(v[2]*a[2,2]+v[1]*a[1,2]+v[0]*a[0,2])+v[1]*(a[1,2]*v[2]+v[1]*A[1,1]+v[0]*a[0,1])+v[0]*(a[0,2]*v[2]+a[0,1]*v[1]+v[0]*a[0,0])
 return (
	  v[2]*(v[2]*a[2][2]+v[1]*a[1][2]+v[0]*a[0][2])
	+ v[1]*(a[1][2]*v[2]+v[1]*a[1][1]+v[0]*a[0][1])
	+ v[0]*(a[0][2]*v[2]+a[0][1]*v[1]+v[0]*a[0][0]) );
}

double FastOps::btAb_Asym_3x3 (const vnl_matrix_fixed<double,3,3>& A, const vnl_vector_fixed<double,3>& b)
{
    double const* a = A.data_block();
    double const* v = b.data_block();

    return (
         v[2]*(v[2]*a[8]+v[1]*a[5]+v[0]*a[2])
       + v[1]*(a[5]*v[2]+v[1]*a[4]+v[0]*a[1])
       + v[0]*(a[2]*v[2]+a[1]*v[1]+v[0]*a[0]) );
}*/

}
