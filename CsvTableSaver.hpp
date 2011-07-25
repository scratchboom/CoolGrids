#pragma once

class CsvTableSaver{
private:
	std::stringstream captionsLine;
	std::stringstream content;
	std::set<std::string> captionsSet;
	int count;
	int currCount;
	bool firstInLine;
	bool captionsAdded;
public:

	CsvTableSaver(){
		count = 0;
		currCount = 0;
		captionsAdded = false;
	}

	void addRow(const std::string& cap,double val){

		if (!captionsAdded) {

			bool found = captionsSet.find(cap)!=captionsSet.end();
			if (!found) {

				if (currCount != 0)
					captionsLine << ",";

				captionsLine << '"' << cap << '"';
				captionsSet.insert(cap);
			}else {
				captionsAdded = true;
				count = captionsSet.size();
			}
		}

		if (captionsAdded && currCount == count) {
			content << std::endl;
			currCount = 0;
		}

		if(currCount!=0)content  << ",";
		content << val;
		currCount++;
	}

#define ADD_ROW(val) addRow(#val,val)

	void save(const std::string& filename) {
		std::ofstream out(filename.c_str());
		out << captionsLine.str() << std::endl;
		out << content.str();
		out.close();
	}

};

