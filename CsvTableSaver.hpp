#pragma once

class CsvTableSaver{
private:
	std::stringstream content;
	bool captionsAdded;
public:

	CsvTableSaver(){
		captionsAdded = false;
	}

	void addRow(const std::string& cap1,double val1){
		if(!captionsAdded){
			content << "\"" << cap1 << "\"";
			captionsAdded = true;
			content << std::endl;
		}

		content << val1 ;
		content << std::endl;
	}

	void addRow(const std::string& cap1, double val1,
			    const std::string& cap2, double val2) {
		if (!captionsAdded) {
			content << "\"" << cap1 <<  "\"";
			content << ",";
			content << "\"" << cap2 <<  "\"";
			content << std::endl;
			captionsAdded = true;
		}

		content << val1;
		content << ",";

		content << val2;

		content << std::endl;
	}

	void addRow(const std::string& cap1, double val1,
			    const std::string& cap2, double val2,
			    const std::string& cap3, double val3) {
		if (!captionsAdded) {
			content << "\"" << cap1 << "\"";
			content << ",";
			content << "\"" << cap2 << "\"";
			content << ",";
			content << "\"" << cap3 << "\"";
			content << std::endl;
			captionsAdded = true;
		}

		content << val1;
		content << ",";
		content << val2;
		content << ",";
		content << val3;

		content << std::endl;
	}

	void addRow(const std::string& cap1, double val1,
			    const std::string& cap2, double val2,
			    const std::string& cap3, double val3,
			    const std::string& cap4, double val4) {
		if (!captionsAdded) {
			content << "\"" << cap1 << "\"";
			content << ",";
			content << "\"" << cap2 << "\"";
			content << ",";
			content << "\"" << cap3 << "\"";
			content << ",";
			content << "\"" << cap4 << "\"";
			content << std::endl;
			captionsAdded = true;
		}

		content << val1;
		content << ",";
		content << val2;
		content << ",";
		content << val3;
		content << ",";
		content << val4;

		content << std::endl;
	}

	void addRow(const std::string& cap1, double val1,
			    const std::string& cap2, double val2,
			    const std::string& cap3, double val3,
			    const std::string& cap4, double val4,
			    const std::string& cap5, double val5) {
		if (!captionsAdded) {
			content << "\"" << cap1 << "\"";
			content << ",";
			content << "\"" << cap2 << "\"";
			content << ",";
			content << "\"" << cap3 << "\"";
			content << ",";
			content << "\"" << cap4 << "\"";
			content << ",";
			content << "\"" << cap5 << "\"";
			content << std::endl;
			captionsAdded = true;
		}

		content << val1;
		content << ",";
		content << val2;
		content << ",";
		content << val3;
		content << ",";
		content << val4;
		content << ",";
		content << val5;

		content << std::endl;
	}

#define ADD_ROW(expr1) addRow(#expr1,expr1)

#define ADD_ROW(expr1,expr2) addRow(#expr1,expr1,\
		                            #expr2,expr2)

#define ADD_ROW(expr1,expr2,expr3) addRow(#expr1,expr1,\
		                                  #expr2,expr2,\
		                                  #expr3,expr3)

#define ADD_ROW(expr1,expr2,expr3,expr4) addRow(#expr1,expr1,\
		                                        #expr2,expr2,\
		                                        #expr3,expr3,\
		                                        #expr4,expr4)

#define ADD_ROW(expr1,expr2,expr3,expr4,expr5) addRow(#expr1,expr1,\
		                                              #expr2,expr2,\
		                                              #expr3,expr3,\
		                                              #expr4,expr4,\
		                                              #expr5,expr5)


	void save(const std::string& filename) {
		std::ofstream out(filename.c_str());
		out << content.str();
		out.close();
	}

};

