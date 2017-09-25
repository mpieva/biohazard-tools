#include <Judy.h>

void finalize_diet( Pvoid_t *pjl )
{
    Word_t i = 0 ;
    Pvoid_t pval = JudyLFirst( *pjl, &i, 0 ) ; 
    while( pval )
    {
        Judy1FreeArray( pval, 0 ) ;
        pval = JudyLNext( *pjl, &i, 0 ) ;
    }
    JudyLFreeArray( pjl, 0 ) ;
}

