/*
 * io.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include <iostream>
#include <ostream>

#include "data_model.h"
#include "global_variables.h"
#include "io.h"

using namespace std;

// I/O:

ostream& operator<< (ostream& os, const link_type& lt) {
    os << et2label[lt.et1]
        << " " << rat2label[lt.rat13] << " "
        << et2label[lt.et3];
    return os;
}

ostream& operator<< (ostream& os, const event_type& evt) {
    os << ec2label[evt.ec]
            << " \"" << et2label[evt.et1]
            << " " << rat2label[evt.rat13] << " "
            << et2label[evt.et3]
            << "\"";
    return os;
}

ostream& operator<< (ostream& os, const event& ev) {
    if (ev.e1 < 0) {
        os << ec2label[ev.ec]
                << " \"some " << et2label[-ev.e1]
                << " " << rat2label[ev.rat13] << " some "
                << et2label[-ev.e3]
                << "\"";
    } else {
        os << ec2label[ev.ec]
                << " \"" << e2label[ev.e1]
                << " " << rat2label[ev.rat13] << " "
                << e2label[ev.e3]
                << "\"";
    }
    return os;
}

void dump_data () // dump ev2data to console for debugging
{
    cout << "!" << endl;
    for (auto& [ev2, evd2] : ev2data)
        cout << ev2 << " " << evd2.attempt_rate << " " << evd2.t << endl;
}


