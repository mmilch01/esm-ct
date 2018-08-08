/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#ifndef _COMMON_INCLUDED_

#define _COMMON_INCLUDED_
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#ifndef _WIN32
#ifndef _4DFP
#define _4DFP
#endif
#endif

#ifdef _4DFP
extern unsigned char    __ctype[];
#endif

#define MAXL 256
#define ARRAY_MAXDIMS 5

#ifndef mmin
#define mmin(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef mmax
#define mmax(a,b) (((a)>(b))?(a):(b))
#endif

#ifdef __GNUG__
#else 
#define NULL 0
#endif //__GNUG__
bool IsBigEndian();
#endif //_COMMON_INCLUDED_
