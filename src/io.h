// make sure this file is only included once:
#ifndef INC_IO_H
#define INC_IO_H

#include "data_model.h"

namespace tricl {

ostream& operator<< (ostream& os, const link_type& lt);

ostream& operator<< (ostream& os, const event_type& evt);

ostream& operator<< (ostream& os, const angle_type& at);

ostream& operator<< (ostream& os, const event& ev);

ostream& operator<< (ostream& os, const event_data& evd);

}

using tricl::operator<<;

void log_state ();

void read_links_csv (
        string filename,
        int skip_rows,
        int max_rows,
        char delimiter,
        int e1_col,
        int rat13_col,
        int e3_col,
        entity_type et1_default,
        relationship_or_action_type rat13_fixed,
        entity_type et3_default,
        string e1_prefix,
        string e3_prefix
        );

void dump_links ();

void dump_data ();

#endif

