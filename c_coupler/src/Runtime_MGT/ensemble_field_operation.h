/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef ENSEMBLE_FIELD_OPERATION
#define ENSEMBLE_FIELD_OPERATION

#include <vector>
#include "memory_mgt.h"


enum {
	ENSEMBLE_OP_TYPE_SUM,
	ENSEMBLE_OP_TYPE_MIN,
	ENSEMBLE_OP_TYPE_MAX,
	ENSEMBLE_OP_TYPE_MEAN,
	ENSEMBLE_OP_TYPE_ANOMALY,
	ENSEMBLE_OP_TYPE_ANY
};


class Member_to_set_operation
{
	private:
		std::vector<Field_mem_info*> member_fields_inst;
		Field_mem_info *set_field_inst;
		void **member_fields_data_buffer;
		void *set_field_data_buffer;
		int operation_type;
		int specified_member_index;
		int field_size;
		char *data_type;                       // 0: float; 1: double; 2: others
		
	public:
		Member_to_set_operation(std::vector<Field_mem_info*>&, Field_mem_info*, int, int);
		~Member_to_set_operation();
		void execute();
};

#endif

