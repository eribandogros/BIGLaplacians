// SparseMatrix.h: interface for the CSparseMatrix class.
// already modify the mulTransMatMat into linear time
//////////////////////////////////////////////////////////////////////

#ifndef _SPARSE_MATRIX_
#define _SPARSE_MATRIX_

#include <stdio.h>
#include <assert.h>
#include <math.h>

#if defined(__linux__) || defined(__APPLE__)

#define fscanf_s fscanf

#endif

#if 0
#define  ___DEBUG_YZ
#endif

//#include "stdafx.h"

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif


// User-defined tolerancy
#define TOL 0.00005
#define	ZERO_TOL 1e-12

class CMatrixElement 
{
public:
	int i;
	int j;
	double value;
	CMatrixElement *rowNext;
	CMatrixElement *colNext;

	CMatrixElement(int newi=0, int newj=0, double newValue=0.0)
	{
		i = newi;
		j = newj;
		value = newValue;
		rowNext = colNext = 0;
	}

	//virtual ~CMatrixElement() { if (rowNext != 0) delete rowNext; }

    //~CMatrixElement() { if (rowNext != 0) delete rowNext; }
};

class CSparseMatrix
{
	//protected:
public:

	int numRows;
	int numCols;
	CMatrixElement* *rowList;
	CMatrixElement* *colList;
	double* diagonal;


public:

	CSparseMatrix(int nRows)
	{
		numRows = numCols = 0;
		rowList = colList = NULL;
		diagonal = NULL;
		setDimensions(nRows);
	}

	CSparseMatrix(int nRows, int nCols) //
	{
		numRows = numCols = 0;
		rowList = colList = NULL;
		diagonal = NULL;
		setDimensions(nRows,nCols);
	}

	CSparseMatrix(int nRows, 
		int numEl, 
		int i[], 
		int j[], 
		double vals[])
	{
		numRows = numCols = 0;
		rowList = colList = NULL;
		diagonal = NULL;
		setDimensions(nRows);
		setValues(numEl,i,j,vals);
	}

	~CSparseMatrix()
	{
		Cleanup();
	}

    void CleanRow(int i)
    {
        CMatrixElement * theElm;
        CMatrixElement * nextElm;
        theElm = rowList[i];
        nextElm = theElm->rowNext;
        for (;nextElm != NULL;)
        {
            if (theElm != NULL)
            {
                delete theElm;
            }
            theElm = nextElm;
            nextElm = nextElm->rowNext;
        }
        if (theElm != NULL)
        {
            delete theElm;
        }
    }

    void Cleanup()
    {
        // Only delete rows since the matrix element deletes rowNext
        for(int i = 0; i < numRows; i++)
            if(rowList[i] != NULL)
                CleanRow(i);
        
        
        if(rowList != NULL)  delete [] rowList;  rowList  = NULL;
        if(colList != NULL)  delete [] colList;  colList  = NULL;
        if(diagonal != NULL) delete [] diagonal; diagonal = NULL;
    }

	void Cleanup_old()
	{
		// Only delete rows since the matrix element deletes rowNext
		for(int i = 0; i < numRows; i++)
			if(rowList[i] != NULL)
				delete rowList[i];

		if(rowList != NULL)  delete [] rowList;  rowList  = NULL;
		if(colList != NULL)  delete [] colList;  colList  = NULL;
		if(diagonal != NULL) delete [] diagonal; diagonal = NULL;
	}

    void CopyFromMatrix(CSparseMatrix * A)
    {
        Cleanup();
        numRows = numCols = 0;
        rowList = colList = NULL;
        diagonal = NULL;
        setDimensions(A->numRows,A->numCols);

        CMatrixElement *theElem;
        for(int i = 0; i < A->numRows; i++)
            for(theElem = A->rowList[i]; theElem != NULL; theElem = theElem->rowNext)
                set1Value(theElem->i, theElem->j, theElem->value);
    }

	void
		 setValues(int numEl, int i[], int j[], double vals[])
	{
		for(int idx = 0; idx < numEl; idx++)
			set1Value(i[idx],j[idx],vals[idx]);
	}

	void
		 set1Value(int i, int j, double val)
	{
		// Insertion in rows
		CMatrixElement *theElem = new CMatrixElement(i,j,val);
		theElem->rowNext = rowList[i];
		rowList[i] = theElem;

		// Insertion in columns
		theElem->colNext = colList[j];
		colList[j] = theElem;

		// If on the diagonal, store it for fast access
		if(i==j)
		{
			diagonal[i] = val;
		}
	}

	void
		 modify1Value(int i, int j, double val)
	{
		CMatrixElement *theElem = GetElement(i,j);

		if(theElem == NULL)
			set1Value(i,j,val);
		else
			theElem->value = val;
	}

	int
		 DeleteElement(int i, int j)
		//not fully tested
	{
        //guess there is memory leak... M(i,j) is not deleted. 

		CMatrixElement *theElem;
		CMatrixElement *leftElem, *rightElem;
		CMatrixElement *aboveElem, *underElem;
		leftElem = rightElem = aboveElem = underElem = NULL;

		for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
		{
			if(theElem->j == j)
			{
				rightElem = theElem->rowNext;
				break;
			}
			leftElem = theElem;
		}

		if (theElem->j != j)
			return -1;
		for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
		{
			if(theElem->i == i)
			{
				underElem = theElem->colNext;
				break;
			}
			aboveElem = theElem;
		}
		if (theElem->i != i)
			return -1;
		if (leftElem == NULL) // Means A(i,j) is the first entry of a row
		{
			rowList[i] = rightElem;
		}
		else // not nessesary to chek rightElem is NULL
		{
			leftElem->rowNext = rightElem;
		}

		if (aboveElem == NULL) // Means A(i,j) is the first entry of a col
		{
			colList[j] = underElem;
		}
		else // not nessesary to chek rightElem is NULL
		{
			aboveElem->colNext = underElem;
		}
		if (i == j)
		{
			diagonal[i] = 0;
		}
		return 1;
	}

    int  GetElementNumber()
        
    {
        int count = 0;
        CMatrixElement *theElem = NULL;
        for(int j = 0; j < numCols; j++)
        {
            for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
                count++;
        }
        return count;
    }

    int  GetNonEmptyRowNum()

    {
        int count = 0;
        for(int j = 0; j < numRows; j++)
        {
            if (rowList[j] != NULL)
            {
                count++;
            }
        }
        return count;
    }

    int  GetMaxColNumer()

    {
        int count = 0;
        int temp = 0;
        CMatrixElement *theElem = NULL;
        for(int j = 0; j < numCols; j++)
        {
            count = 0;
            for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
                count++;
            if ( count > temp)
                temp = count;
        }
        return temp;
    }


	void add1Value(int i, int j, double val)
	{
		CMatrixElement *theElem = GetElement(i,j);

		if(theElem == NULL)
		{
			if ( fabs(val) > ZERO_TOL)
			{
				set1Value(i,j,val);
			}
		}
		else
		{
			theElem->value += val;
			if ( fabs(theElem->value) < ZERO_TOL)
			{
				DeleteElement(i,j);
			}
		}

	}

	void addOneValue(int i, int j, double val)
	{
		CMatrixElement *theElem = GetElement(i,j);

		if(theElem == NULL)
			set1Value(i,j,val);
		else
			theElem->value += val;
	}

	void setRow(int i, CMatrixElement *head)
	{
		// Set it in the row
		rowList[i] = head;
		// And in the column (and diagonal)
		CMatrixElement *theElem;
		for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
		{
			theElem->colNext = colList[theElem->j];
			colList[theElem->j] = theElem;
			if(i == theElem->j)
				diagonal[i] = theElem->value;
		}
	}

	void setDimensions(int nRows)
	{
		// Clean up anyway. Safer, since life is a jungle.
		Cleanup();
		numRows = nRows;
		numCols = nRows;
		/*
		if(rowList)
		delete [] rowList;
		if(colList)
		delete [] colList;
		if(diagonal)
		delete [] diagonal;*/
		rowList = new CMatrixElement*[numRows];
		colList = new CMatrixElement*[numCols];
		diagonal = new double[numRows];
		for(int k = 0; k < numRows; k++)
		{
			diagonal[k] = 0.;
			rowList[k] = colList[k] = NULL;
		}	

	}

	void setDimensions(int nRows, int nCols)
	{
		// Clean up anyway. Safer, since life is a jungle.
		Cleanup();
		numRows = nRows;
		numCols = nCols;
		rowList = new CMatrixElement*[numRows];
		colList = new CMatrixElement*[numCols];
		diagonal = new double[numRows];
		for(int k = 0; k < numRows; k++)
		{
			diagonal[k] = 0.;
			rowList[k] = NULL;
		}
		for(int l = 0; l < numCols; l++)
			colList[l] = NULL;
	}

	CMatrixElement*
		 GetElement(int i, int j)
	{
		CMatrixElement *theElem;
		for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
			if(theElem->j == j)
				return theElem;
		return NULL;
	}

	double
		 GetValue(int i, int j)
	{
		CMatrixElement *theElem;
		for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
			if(theElem->j == j)
				return theElem->value;
		return 0.0;
	}

	double
		 diagonalElement(int i)
	{
		assert(i < numRows);
		return diagonal[i];
	}

	void
		 Print()
	{
		CMatrixElement *theElem;
		for(int i = 0; i < numRows; i++)
			for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
				printf("i=%d, j=%d: %f\n", theElem->i, theElem->j, theElem->value);
	}

    void
        PrintDense()
    {
        // print column index
        printf("idx");
        for (int j = 0; j < numCols; j++)
        {
            printf("\t|   %d", j);
        }
        printf("\t|\n");
        for(int i = 0; i < numRows; i++)
        {
            // print row index
            printf("%d", i);
            for (int j = 0; j < numCols; j++)
            {
                if (GetElement(i,j) == NULL)
                {
                    printf("\t|");
                }
                else
                {
                    printf("\t|   %.0f", GetValue(i,j));
                }
            }
            printf("\t|\n");
        }
    }

	void
		 PrintMathematica(FILE *fp)
	{
		int i, j;
		fprintf(fp,"m = {");
		for(i = 0; i < numRows; i++)
		{
			fprintf(fp,"\n{");
			for(j = 0; j < numCols; j++)
			{
				fprintf(fp,"%f",GetValue(i,j));
				if(j != numCols-1) fprintf(fp,", ");
			}
			fprintf(fp,"}");
			if(i != numRows-1) fprintf(fp,",");
		}
		fprintf(fp,"}\n\n");
	}
	void 
		 PrintMathematica_wyz(FILE *fp)
	{
		int i, j;
		fprintf(fp,"This is a %d by %d matrix.\n",numRows,numCols);
		fprintf(fp,"m = {");
		for(i = 0; i < numRows; i++)
		{
			fprintf(fp,"\n{");
			for(j = 0; j < numCols; j++)
			{
				double a = GetValue(i,j);
				if( GetElement(i,j) != NULL)
				{
					fprintf(fp,"%d,%d %e",i,j,a);
					if(j != numCols-1) fprintf(fp,", ");
				}
			}
			fprintf(fp,"}");
			if(i != numRows-1) fprintf(fp,",");
		}
		fprintf(fp,"}\n\n");
	}

	void 
		 PrintMathematica_wyz2(FILE *fp)
	{
        int i, j;
        fprintf(fp,"This is a %d by %d matrix.\n",numRows,numCols);
        fprintf(fp,"m = {");
        for(i = 0; i < numRows; i++)
        {
            fprintf(fp,"\n{");
            for(j = 0; j < numCols; j++)
            {
                double a = GetValue(i,j);
                if ( fabs(a) < 0.0001)
                    fprintf(fp,"  ");
                else
                {
                    if (a > 0)
                        fprintf(fp,"+%d",1);
                    else
                        fprintf(fp,"-%d",1);
                }
            }
            fprintf(fp,"}");
            if(i != numRows-1) fprintf(fp,",");
        }
        fprintf(fp,"}\n\n");
	}
	void
		PrintVectorMathematica(FILE *fp, double theVec[], int n)
	{
		int i;
		fprintf(fp,"v = {");
		for(i = 0; i < n; i++)
		{
			fprintf(fp,"%f",theVec[i]);
			if(i != n-1) fprintf(fp,", ");
		}
		fprintf(fp,"}\n\n");
	}


	void
		 multMatVec(double *src,
		double *dest)
	{
		assert(src && dest);
		CMatrixElement *theElem = NULL;
		for(int i = 0; i < numRows; i++)
		{
			double sum = 0;
			for(theElem = rowList[i];
				theElem != NULL;
				theElem = theElem->rowNext)
				sum += theElem->value * src[theElem->j];
			dest[i] = sum;
		}
	}
	void 	 multMatVec_yz(double *src,
		double * &det)
	{
		double *dest;
		dest = new double [numRows];
		CMatrixElement *theElem = NULL;
		for(int i = 0; i < numRows; i++)
		{
			double sum = 0;
			for(theElem = rowList[i];
				theElem != NULL;
				theElem = theElem->rowNext)
				sum += theElem->value * src[theElem->j];
			dest[i] = sum;
		}
		//for(int i = 0; i < numRows; i++)
		//{
		//	det[i] = dest[i];
		//}
		det = dest;
		//delete [] dest;
	}

	void
		 multTransMatVec(double *src,
		double *dest)
	{
		assert(src && dest);
		double sum;

		CMatrixElement *theElem = NULL;
		for(int j = 0; j < numCols; j++)
		{
			sum = 0.0;
			for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
				sum += theElem->value * src[theElem->i];
			dest[j] = sum;
		}
	}

	void
		 multTransMatVec_yz(double *src,
		double *det)
	{
		
		double sum;
		double * dest;
		dest = new double[numCols];
		CMatrixElement *theElem = NULL;
		for(int j = 0; j < numCols; j++)
		{
			sum = 0.0;
			for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
				sum += theElem->value * src[theElem->i];
			dest[j] = sum;
			det[j] = dest[j];
		}
		delete [] dest;
	}

	//void
	// Transpose()
	//{
	//
	//	CSparseMatrix* tempMat = new CSparseMatrix(numCols,numRows);
	//
	//}

	void
		 multTransMatMat_yz()
	{
		// M = transpose(M)*M  <-> A * B
		CSparseMatrix* tempMat = new CSparseMatrix(numCols);

		int j;
		CMatrixElement *theElem;
		
		#ifdef ___DEBUG_YZ
			FILE * fp;
			fp = fopen("mat_mup_log.txt","w");
			fclose(fp);
		#endif
		
		for(j = 0; j < numCols; j++)
		{
			theElem = colList[j];
			if (theElem == NULL)
				continue;

			for(CMatrixElement * theElem_e = colList[j]; theElem_e != NULL; theElem_e = theElem_e->colNext)  //D^T: col means A.row
			{
				int r_id = theElem_e->i;
				int tem_rid = theElem_e->j;
				for(CMatrixElement *theElem_r = rowList[r_id]; theElem_r != NULL; theElem_r = theElem_r->rowNext)
				{
					#ifdef ___DEBUG_YZ
						fp = fopen("mat_mup_log.txt","a+");
					#endif

					int tem_cid = theElem_r->j;

					#ifdef ___DEBUG_YZ
						fprintf(fp,"j=%d r_id=%d tem_rid=%d tem_cid=%d\n",j,r_id,tem_rid,tem_cid);
					#endif

					double value = theElem_e->value * theElem_r->value;

					#ifdef ___DEBUG_YZ
						fprintf(fp,"A[%d,%d] = %f, B[%d,%d] = %f, C[%d,%d] += %f\n",tem_rid,r_id,theElem_e->value,r_id,tem_cid,theElem_r->value,tem_rid,tem_cid,value);
					#endif

					tempMat->add1Value(tem_rid,tem_cid,value);

					#ifdef ___DEBUG_YZ
						fclose(fp);
					#endif

				}
			}

		}
		// copy tempMat to this
		//Cleanup();
		setDimensions(numCols);
		rowList = tempMat->rowList;
		colList = tempMat->colList;
		diagonal = tempMat->diagonal;
		// delete tempMat; // must fix that ! -> memory leaks here
	}

	void
		 writeToFile(FILE *fp)
	{
		CMatrixElement *theElem;
		for(int i = 0; i < numRows; i++)
			for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
				fprintf(fp, "%d %d %lf\n", theElem->i, theElem->j, theElem->value);
	}

	void
		 readFromFile(FILE *fp)
	{
		int i,j;
		double value;
		while (fscanf_s(fp, "%d %d %lf", &i,&j,&value) != EOF)
		{
			set1Value(i,j,value);
		}

	}



	void
		 multTransMatMat()
	{
		// M = transpose(M)*M
		CSparseMatrix* tempMat = new CSparseMatrix(numCols);
		double *colVec = new double[numRows];
		double *des_colVec = new double[numCols];

		int i, j;
		CMatrixElement *theElem;

		for(j = 0; j < numCols; j++)
		{
			theElem = colList[j];
			if (theElem == NULL)
				continue;
			// initialize the result
			for(i = 0; i < numCols; i++)
				des_colVec[i] = 0.0;

			// grab column j
			for(i = 0; i < numRows; i++)
				colVec[i] = 0.0;
			for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
				colVec[theElem->i] = theElem->value;

			// compute a column of M
			multTransMatVec(colVec, des_colVec);

			// store in tempMat
			for(i = 0; i < numCols; i++)
				if(fabs(des_colVec[i]) > 0.0001)
					tempMat->set1Value(i,j,des_colVec[i]);
		}

		// copy tempMat to this
		// Cleanup();
		setDimensions(numCols);
		rowList = tempMat->rowList;
		colList = tempMat->colList;
		diagonal = tempMat->diagonal;
		delete [] colVec;
		delete [] des_colVec;
		// delete tempMat; // must fix that ! -> memory leaks here
	}

	void
		 AddMatrix(CSparseMatrix *mat)
	{
		int i;
		CMatrixElement *theElem, *matElem;
		for(i = 0; i < numRows; i++)
		{
			for(matElem = mat->rowList[i]; matElem != NULL; matElem = matElem->rowNext)
			{
				theElem = GetElement(matElem->i,matElem->j);
				if(theElem == NULL)
				{
					theElem = new CMatrixElement(matElem->i,matElem->j,matElem->value);
					theElem->rowNext = rowList[i];
					rowList[i] = theElem;
					theElem->colNext = colList[theElem->j];
					colList[theElem->j] = theElem;
				}
				else
					theElem->value += matElem->value;
			}
		}
		for(i = 0; i < numRows; i++)
			diagonal[i] += mat->diagonal[i];
	}

	void
		 ScaleRow(int i, double s)
	{
		CMatrixElement *theElem;
		for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
			theElem->value *= s;
		diagonal[i] *= s;
	}

    //matrix multiplication: result = this * mat
    CSparseMatrix *
        MultMatrix_bb(CSparseMatrix *mat){
            //check if the size of matrices corresponds
            if(this->numCols != mat->numRows)
                return NULL;

            CSparseMatrix *result =new CSparseMatrix(numRows,mat->numCols);
            CMatrixElement *theElem,*matElem;
            for(int i=0; i<this->numRows; i++){
                for(theElem = this->rowList[i]; theElem != NULL; theElem = theElem->rowNext){
                    int k = theElem->j;
                    for(matElem = mat->rowList[k]; matElem != NULL; matElem = matElem->rowNext){
                        int j=matElem->j;
                        result->add1Value(i,j,theElem->value * matElem->value);

                    }
                }

            }

            for(int i=0; i<result->numRows; i++)
                result->diagonal[i] = result->GetValue(i,i);

            return result;

            ////test
            //TRACE("left matrix:\n");
            //this->Trace();
            //TRACE("\n");
            //TRACE("right matrix:\n");
            //mat->Trace();
            //TRACE("result matrix:\n");
            //result->Trace();



    }

    CSparseMatrix * 
        Scale_Matrix(double scale){
            CSparseMatrix * result = new CSparseMatrix(numRows,numCols);
            CMatrixElement *theElem;
            for(int i=0;i<numRows;i++){
                for(theElem = rowList[i];theElem!=NULL;theElem = theElem->rowNext)
                    result->add1Value(i,theElem->j,scale*(theElem->value));
            }

            for(int i=0; i<result->numRows; i++)
                result->diagonal[i] =result->GetValue(i,i); 

            return result;
    }

    void
        MinusMatrix(CSparseMatrix *mat)
    {
        int i;
        CMatrixElement *theElem, *matElem;
        for(i = 0; i < numRows; i++)
        {
            for(matElem = mat->rowList[i]; matElem != NULL; matElem = matElem->rowNext)
            {
                theElem = GetElement(matElem->i,matElem->j);
                if(theElem == NULL)
                {
                    theElem = new CMatrixElement(matElem->i,matElem->j,-1*(matElem->value));
                    theElem->rowNext = rowList[i];
                    rowList[i] = theElem;
                    theElem->colNext = colList[theElem->j];
                    colList[theElem->j] = theElem;
                }
                else
                    theElem->value -= matElem->value;
            }
        }
        for(i = 0; i < numRows; i++)
            diagonal[i] -= mat->diagonal[i];
    }


    CSparseMatrix *
        SetTranspose(){

            CSparseMatrix *matrix_tran = new CSparseMatrix(numCols,numRows);
            CMatrixElement *theElem;
            for(int i=0; i<numRows; i++){
                for(theElem = rowList[i]; theElem !=NULL; theElem=theElem->rowNext){
                    int j = theElem->j;
                    matrix_tran->add1Value(j,i,theElem->value);
                }
            }
            return matrix_tran;

    }


	//***************************************
	// preconditionedBiConjugateGradient
	//***************************************
	unsigned int 
		preconditionedBiConjugateGradient(double x[],
		double b[],
		double tol,
		const unsigned int iter_max)
	{
		double *dr = new double[numRows];
		double *drb = new double[numRows];
		double *dp = new double[numRows];
		double *dpb = new double[numRows];
		double *dz = new double[numRows];
		double *dzb = new double[numRows];
		double *dAp = new double[numRows];
		double *dATpb = new double[numRows];

		assert(dr && drb && dp && dpb && dz && dzb && dAp && dATpb);
		double mag_r, mag_rOld, mag_pbAp, mag_Residual, Residual0, alpha, beta;

		multMatVec(x,dAp);
		mag_r = mag_Residual = Residual0 = 0.;
		int i = 0;
		for(i = 0; i < numRows; i++)
		{

			dr[i] = drb[i] = b[i] - dAp[i];
			dp[i] = dpb[i] = dz[i] = dzb[i] = dr[i]/diagonalElement(i);	// Simple preconditioning
			mag_r += drb[i] * dz[i];
			mag_Residual += dz[i] * dz[i];
			Residual0 += b[i]*b[i]/(diagonalElement(i)*diagonalElement(i));
		}

		mag_Residual = Residual0*100; // Force the first iteration anyway.
		if(Residual0 == 0)
			Residual0 = 1.;	// To make it work even if ||b|| = 0
		unsigned int nbIter = 0;
		while(mag_Residual > tol && nbIter < iter_max)
		{
			nbIter++;
			multMatVec(dp,dAp);
			multTransMatVec(dpb,dATpb);
			mag_pbAp = 0.0;
			for(i = 0; i < numRows; i++)
				mag_pbAp += dpb[i] * dAp[i];

			if(mag_pbAp == 0)
			{
				//fprintf(stderr,"OOOOOCH!!! (mag_pbAp==0)\n");
				//return 0;
			}

			if(mag_r == 0 && mag_pbAp == 0)
				alpha = 1;
			else
				alpha = mag_r / mag_pbAp;
			mag_rOld = mag_r;
			mag_r = 0.0;
			for(i = 0; i < numRows; i++)
			{
				x[i] += alpha * dp[i];
				dr[i] -= alpha * dAp[i];
				drb[i] -= alpha * dATpb[i];
				dz[i] = dr[i]/diagonalElement(i);
				dzb[i] = drb[i]/diagonalElement(i);
				mag_r += drb[i] * dz[i];
			}

			if(mag_rOld == 0)
			{
				//fprintf(stderr,"OOOOOCH!!! (mag_rOld==0)\n");
				//return 0;
			}

			if(mag_r == 0 && mag_rOld == 0)
				beta = 1.0;
			else
				beta = mag_r / mag_rOld;
			mag_Residual = 0.;
			for(i = 0; i < numRows; i++)
			{
				dp[i] = dz[i] + beta * dp[i];
				dpb[i] = dzb[i] + beta * dpb[i];
				mag_Residual += dz[i] * dz[i];
			}
		}
		delete [] dr;
		delete [] drb;
		delete [] dp;
		delete [] dpb;
		delete [] dz;
		delete [] dzb;
		delete [] dAp;
		delete [] dATpb;

		return nbIter;
	}


	void
		 Debug()
	{
		int i, j;
		fprintf(stderr,"M [ %i %i ] =\n",numRows,numCols);
		for(i = 0; i < numRows; i++)
		{
			for(j = 0; j < numCols; j++)
				fprintf(stderr,"%f ",GetValue(i,j));
			fprintf(stderr,"\n");
		}
	}

};

#endif // _SPARSE_MATRIX_

