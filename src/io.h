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

void log_status ();

void read_links_csv (
        string filename,    ///< name of infile
        int skip_rows,      ///< no. of rows to skip at start of file (e.g. header lines)
        int max_rows,       ///< no. of rows to read at most
        char delimiter,     ///< delimiter, e.g. " " or "," or "\t"
        int source_col,     ///< no. of column containing e1 labels
        int type_col,       ///< no. of column containing rat13 labels (or -1 if type is fixed to value of "rat")
        int target_col,     ///< no. of column containing e3 labels
        entity_type et1_default, ///< default entity-type for previously unregistered e1s
        relationship_or_action_type fixed_type, ///< fixed type if type_col == -1 (otherwise -1)
        entity_type et3_default, ///< default entity-type for previously unregistered e3s
        string label_prefix ///< prefix to prepend to entity labels before storing them
        );

void dump_links ();

void dump_data ();

#endif

