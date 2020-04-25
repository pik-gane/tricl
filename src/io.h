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
        string e1_prefix,
        int et1_col,
        entity_type et1_default,

        int rat13_col,
        relationship_or_action_type rat13_fixed,

        int e3_col,
        string e3_prefix,
        int et3_col,
        entity_type et3_default
        );

void dump_links ();

void dump_data ();

#endif

