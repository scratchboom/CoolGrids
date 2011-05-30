#pragma once

std::string trim(const std::string& str)
{
  std::string::size_type pos1 = str.find_first_not_of(" \t\r\n");
  std::string::size_type pos2 = str.find_last_not_of(" \t\r\n");
  return str.substr(pos1 == std::string::npos ? 0 : pos1,pos2 == std::string::npos ? str.length() - 1 : pos2 - pos1 + 1);
}

std::string frame(const std::string& filename,double it,const std::string& ext){
	return filename+toZPadString((int)it,5)+"."+ext;
}
