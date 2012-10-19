// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <QtCore>
#include <iostream>
#include <stdio.h>

#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_matrix.h>



extern "C" {

SEXP readLabelFile(SEXP name)
{
	const char *path = CHAR(STRING_ELT(name, 0));
	QFile file(path);
	file.open(QIODevice::ReadOnly);
	
	//std::cout << "reading file " << path << "...\n";
	Rprintf("reading file %s...\n", path);
	QDataStream in(&file);
	
	qint32 dummy, size;
	in >> dummy;
	in >> size;
	
	//std::cout << "file has " << size << " labels\n";
	Rprintf("file has %d labels\n", (int)size);	

	quint8 elem;
	int display;
	SEXP res;
	PROTECT(res=allocVector(INTSXP, size));
	
	for(int i=0; i<size; i++) {
		in >> elem;
		display = elem;
		INTEGER(res)[i] = display;
	}

	UNPROTECT(1);
	return(res);
}

/* grey pixmap object :
$size (nrow and ncol)
$cellres (2 times 1)
$bbox : 0 0 ncol nrow
$bbcent FALSE
$channels "grey"
$grey : matrix with 0s (black) to 1s (white)
!!! we get 255 - 0s from provided data files ! post processed
*/

SEXP readPixmapFile(SEXP name)
{
	const char *path = CHAR(STRING_ELT(name, 0));
	QFile file(path);
	file.open(QIODevice::ReadOnly);
	
	//std::cout << "reading file " << path << "...\n";
	Rprintf("reading file %s...\n", path);
	QDataStream in(&file);

	// get size info
	qint32 dummy, nim, nrow, ncol;
	int inim, inrow, incol;
	in >> dummy >> nim >> nrow >> ncol ;
	inim=nim;
	inrow=nrow;
	incol=ncol;
	//std::cout << "file contains " << nim << " images format " << nrow << "x" << ncol << "\n";
	Rprintf("file contains %d images format %dx%d\n", (int)nim, (int)nrow, (int)ncol);

	// images shall be loaded in a list of matrices, returned raw. conversion to pixmap done externally.
	// first create appropriate R data structure
	SEXP pack, curmat;
	quint8 read;
	int iread;
	PROTECT(pack=allocVector(VECSXP, inim));
	for(int i=0; i<inim; i++) {
		SET_VECTOR_ELT(pack, i, curmat=allocMatrix(INTSXP, inrow, incol));
		for(int j=0; j<inrow; j++) {
			for(int k=0; k<incol; k++) {
				in >> read;
				iread=read;
				INTEGER(curmat)[inrow * k + j] = iread;
			}
		}
	}
	
	UNPROTECT(1);
	return(pack);
}

/*
SEXP getPathsForFormat(SEXP format, SEXP pathIsAbsolute) {
	int pathType = INTEGER(pathIsAbsolute)[0];
     QDir dir("/Users/Pierrick/These/256mod");
     dir.setFilter(QDir::Dirs | QDir::NoDotAndDotDot);
     dir.setSorting(QDir::Name);

     QFileInfoList list = dir.entryInfoList();
     QStringList filters;
     filters << CHAR(STRING_ELT(format, 0));
     //std::cout << "     Bytes Filename" << std::endl;
     //std::cout << "list size is " << list.size() << std::endl;
     
     SEXP nameslist, curnames;
     PROTECT(nameslist = allocVector(VECSXP, list.size()));
	 int labcounts[list.size()];
	 int totsize;
     
     for (int i = 0; i < list.size(); i++) {
         //QFileInfo fileInfo = list.at(i);
         //std::cout << qPrintable(QString("%1").arg(fileInfo.fileName()));
         //std::cout << std::endl;
         QDir curdir = dir;
		QString curpath;
         curdir.cd(list.at(i).fileName());
         curdir.setFilter(QDir::Files);
         curdir.setNameFilters(filters);
		 QFileInfoList files = curdir.entryInfoList();
		 labcounts[i] = files.size();
		 totsize += files.size();
		 SET_VECTOR_ELT(nameslist, i, allocVector(STRSXP, files.size()));
		 curnames = VECTOR_ELT(nameslist, i);
		 for(int j=0; j<files.size(); j++) {
			if(pathType) {
				curpath = files.at(j).absoluteFilePath();
			} else {
				curpath = list.at(i).fileName() + "/" + files.at(j).fileName();
			}
		 	SET_STRING_ELT(curnames, j, mkChar(qPrintable(QString("%1").arg(curpath))));
		 	//std::cout << "processing " << qPrintable(QString("%1").arg(files.at(j).absoluteFilePath())) << std::endl;
		 }

     }

     // build plain vectors
     SEXP outnames, outlab;
     PROTECT(outnames = allocVector(STRSXP, totsize));
     PROTECT(outlab = allocVector(INTSXP, totsize));
     int curind = 0;
     for(int i=0; i<list.size(); i++) {
     	for(int j=0; j<labcounts[i]; j++) {
     		SET_STRING_ELT(outnames, curind, STRING_ELT(VECTOR_ELT(nameslist, i), j));
     		INTEGER(outlab)[curind] = i+1;
     		curind++;
     	}
     }
     
     SEXP out;
     PROTECT(out = allocVector(VECSXP, 2));
     SET_VECTOR_ELT(out, 0, outnames);
     SET_VECTOR_ELT(out, 1, outlab);
     
     UNPROTECT(4);
     return out;
}
*/

/*
void writeMatrix(gsl_matrix *mat) {
	QFile file("/Users/Pierrick/R/log.txt");
	file.open(QIODevice::Append);
	
	QDataStream out(&file);
	// write line to identify output
	out << "logged at" << QDateTime::currentDateTime() << "\n";
	out << "-------------------\n";
	
	int i,j;
	for(i=0; i<mat->size1; i++) {
		for(j=0; j<mat->size2; j++) {
			out << gsl_matrix_get(mat, i, j) << " ";
		}
		out << "\n";
	}
	out << "\n";
}

*/
}
