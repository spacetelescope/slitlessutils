#include <stdio.h>
#include <Python.h>
#include <numpy/npy_no_deprecated_api.h>
#include <numpy/arrayobject.h>
#include "polyclip.h"

static PyObject *_multi(PyObject *self,PyObject *args){

  /* Create objects from the inputs */
  PyObject *lobj,*robj,*bobj,*tobj,*pxobj,*pyobj;
  PyObject *n_polyobj,*poly_indsobj,*xxobj,*yyobj,*nclip_polyobj,*areasobj;
  if (!PyArg_ParseTuple(args, "OOOOOOOOOOOO",&lobj,&robj,&bobj,&tobj,&pxobj,&pyobj,&n_polyobj,&poly_indsobj,&xxobj,&yyobj,&nclip_polyobj,&areasobj)){
    return NULL;
  } 

  
  /* if arrays, then extract them to objects */
  PyObject *larr=PyArray_FROM_OTF(lobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *rarr=PyArray_FROM_OTF(robj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *barr=PyArray_FROM_OTF(bobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *tarr=PyArray_FROM_OTF(tobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *pxarr=PyArray_FROM_OTF(pxobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);
  PyObject *pyarr=PyArray_FROM_OTF(pyobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);
  PyObject *n_polyarr=PyArray_FROM_OTF(n_polyobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *poly_indsarr=PyArray_FROM_OTF(poly_indsobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *xxarr=PyArray_FROM_OTF(xxobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *yyarr=PyArray_FROM_OTF(yyobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *nclip_polyarr=PyArray_FROM_OTF(nclip_polyobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *areasarr=PyArray_FROM_OTF(areasobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);

  /* extract the array data to a C variable */
  int *l = (int*)PyArray_DATA((PyArrayObject*)larr);
  int *r = (int*)PyArray_DATA((PyArrayObject*)rarr);
  int *b = (int*)PyArray_DATA((PyArrayObject*)barr);
  int *t = (int*)PyArray_DATA((PyArrayObject*)tarr);
  float *px = (float*)PyArray_DATA((PyArrayObject*)pxarr);
  float *py = (float*)PyArray_DATA((PyArrayObject*)pyarr);
  int *n_poly=(int*)PyArray_DATA((PyArrayObject*)n_polyarr);
  int *poly_inds=(int*)PyArray_DATA((PyArrayObject*)poly_indsarr);
  int *xx = (int*)PyArray_DATA((PyArrayObject*)xxarr);
  int *yy = (int*)PyArray_DATA((PyArrayObject*)yyarr);
  int *nclip_poly=(int*)PyArray_DATA((PyArrayObject*)nclip_polyarr);
  float *areas = (float*)PyArray_DATA((PyArrayObject*)areasarr);
  
  /* clean up memory */
  Py_DECREF(larr);
  Py_DECREF(rarr);
  Py_DECREF(barr);
  Py_DECREF(tarr);
  Py_DECREF(pxarr);
  Py_DECREF(pyarr);
  Py_DECREF(n_polyarr);
  Py_DECREF(poly_indsarr);
  Py_DECREF(xxarr);
  Py_DECREF(yyarr);
  Py_DECREF(nclip_polyarr);
  Py_DECREF(areasarr);
  
  /* call function */
  int n=n_poly[0];
  polyclip_multi(l,r,b,t,px,py,n,poly_inds,xx,yy,nclip_poly,areas);
  
  //printf("C: %f %f\n",px[0],xx[0]);
  

  n_poly[0]=n;
  
  /* Do something interesting here. */
  Py_RETURN_NONE;
}


static PyObject *_single(PyObject *self,PyObject *args){
  /* Create objects from the inputs */
  PyObject *lobj,*robj,*bobj,*tobj,*nvertsobj;
  PyObject *pxobj,*pyobj,*px_outobj,*py_outobj,*areasobj;
  PyObject *indsobj,*nclip_polyobj,*ri_outobj;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",&lobj,&robj,&bobj,&tobj,&pxobj,&pyobj,&nvertsobj,&px_outobj,&py_outobj,&indsobj,&nclip_polyobj,&areasobj,&ri_outobj)){
    return NULL;
  }

  
  
  /* if arrays, then extract them to objects */
  PyObject *larr=PyArray_FROM_OTF(lobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *rarr=PyArray_FROM_OTF(robj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *barr=PyArray_FROM_OTF(bobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *tarr=PyArray_FROM_OTF(tobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *nvertsarr=PyArray_FROM_OTF(nvertsobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *pxarr=PyArray_FROM_OTF(pxobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);
  PyObject *pyarr=PyArray_FROM_OTF(pyobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);
  PyObject *px_outarr=PyArray_FROM_OTF(px_outobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);
  PyObject *py_outarr=PyArray_FROM_OTF(py_outobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);
  PyObject *areasarr=PyArray_FROM_OTF(areasobj,NPY_FLOAT32,NPY_ARRAY_IN_ARRAY);
  PyObject *indsarr=PyArray_FROM_OTF(indsobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *nclip_polyarr=PyArray_FROM_OTF(nclip_polyobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);
  PyObject *ri_outarr=PyArray_FROM_OTF(ri_outobj,NPY_INT32,NPY_ARRAY_IN_ARRAY);


  
  /* extract the array data to a C variable */
  int *l = (int*)PyArray_DATA((PyArrayObject*)larr);
  int *r = (int*)PyArray_DATA((PyArrayObject*)rarr);
  int *b = (int*)PyArray_DATA((PyArrayObject*)barr);
  int *t = (int*)PyArray_DATA((PyArrayObject*)tarr);
  int *nverts = (int*)PyArray_DATA((PyArrayObject*)nvertsarr);
  float *px = (float*)PyArray_DATA((PyArrayObject*)pxarr);
  float *py = (float*)PyArray_DATA((PyArrayObject*)pyarr);
  float *px_out = (float*)PyArray_DATA((PyArrayObject*)px_outarr);
  float *py_out = (float*)PyArray_DATA((PyArrayObject*)py_outarr);
  float *areas=(float*)PyArray_DATA((PyArrayObject*)areasarr);
  int *inds=(int*)PyArray_DATA((PyArrayObject*)indsarr);
  int *nclip_poly=(int*)PyArray_DATA((PyArrayObject*)nclip_polyarr);
  int *ri_out = (int*)PyArray_DATA((PyArrayObject*)ri_outarr);

  
  /* clean up memory */
  Py_DECREF(larr);
  Py_DECREF(rarr);
  Py_DECREF(tarr);
  Py_DECREF(barr);
  Py_DECREF(nvertsarr);	
  Py_DECREF(pxarr);
  Py_DECREF(pyarr);
  Py_DECREF(px_outarr);
  Py_DECREF(py_outarr);
  Py_DECREF(areasarr);
  Py_DECREF(indsarr);
  Py_DECREF(nclip_polyarr);
  Py_DECREF(ri_outarr);


  //printf("%i %i %i %i %i\n",l[0],r[0],t[0],b[0],nverts[0]);
  
  /* call function */
  int n=nclip_poly[0];
  polyclip_single(l[0],r[0],b[0],t[0],px,py,nverts[0],inds,&n,areas,px_out,py_out,ri_out);
  //printf("%i\n",inds[0]);
  nclip_poly[0]=n;
  
  /* Do something interesting here. */
  Py_RETURN_NONE;
}




static PyMethodDef module_methods[]={
  { "multi", (PyCFunction)_multi, METH_NOARGS,NULL },
  { "multi", _multi, METH_VARARGS, "A python driver to call polyclip_multi.\nA function written by J.D. Smith\n"},
  { "single", (PyCFunction)_single, METH_NOARGS,NULL },
  { "single", _single, METH_VARARGS, "A python driver to call polyclip_single.\nA function written by J.D. Smith\n"},
  { NULL, NULL, 0, NULL }
};




  
static struct PyModuleDef _CPolyClip =
{
    PyModuleDef_HEAD_INIT,
    "cpolyclip", /* name of module */
    NULL,
    -1, 
    module_methods
};

PyMODINIT_FUNC PyInit_cpolyclip(void)
{
  import_array();
  return PyModule_Create(&_CPolyClip);
}








