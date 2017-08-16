#! /usr/bin/env python
import sys
src_dir = sys.argv[1]

ini = open("default.ini","r")

lines = ini.readlines()

ini.close()

f90sbr = open(src_dir+"/filINIReader.f90","w")
f90sbr2 = open(src_dir+"/filCMDReader.f90","w")

headerString = (
	"\n\tuse stateMod"
	+"\n\tuse stateIndex"
	+"\n\tuse variables1"
	+"\n\tuse initCond"
	+"\n\tuse boundMod"
	+"\n\tuse parametersBackground"
	+"\n\tuse filConstants"
	+"\n\tuse errorMod"
	+"\n\tuse potential"
	+"\n\tuse constants"
	+"\n\tuse wdir"
	+"\n\tuse dnams"
	+"\n\tuse fricoptions"
	+"\n\timplicit none")

f90sbr.write("subroutine readini(name,iflag)"
	+headerString
	+"\n\tcharacter(len=100) :: name"
	+"\n\tinteger :: iflag,ifound"
	+"\n\tcharacter(len=100) :: iniStr\n\n\t"
	+"if (iflag==1) then\n\t"
	+"call setIniFilename(ininam)\n\n\t"
	+"else\n\t"
	+"call setIniFilename(trim(name))\n\t"
	+"endif\n\n\t")

f90sbr2.write("subroutine readcmdinput(name,str_len_mid)"
	+headerString
	+"\ninteger :: istr\n"
	+"integer :: str_len_mid\n"
	+"integer :: str_len_max=100\n"
	+"character(len=100) :: name\n"
	+"do istr = 2,str_len_mid\n"
	+"select case(name(1:istr))\n"
	)

conditional_pre = "if(ifound==1 .or. iflag==1)  "

section = ''
for line in lines:
	line = line.strip()
	l = len(line)
	if l > 2:
		if line[0] == '[' and line[l-1] == ']':
			section = line[1:l-1]
			word = "!********** section "+section+" *************\n\n\t"
			f90sbr.write(word)
		else:
			vali = 0
			for i in range(1,l-1):
				if line[i] == '=':
					vali = i+1
					break
			if line[0:2] == 'b_' and vali>0:
				word = ("!**logical**\n\t"
				+"call getValue(ifound,iflag,\'"+section+"\',\'"+line[:vali-1].strip()+"\',"
				+"iniStr,\'"+line[vali:].strip()+"\')\n\t"
				+conditional_pre+"call intStrToLStr(iniStr,iniStr)\n\t"
				+conditional_pre+"read(iniStr,*) "+line[2:vali-1].strip()+"\n\n\t")
				f90sbr.write(word)
				word = ("!**logical**\n"
				+"case(\'"+line[:vali-1].strip()+"=\',\'"+line[2:vali-1].strip()+"=\')\n"
				+"call intStrToLStr(name(istr+1:str_len_max),name(istr+1:str_len_max))\n"
				+"read(name(istr+1:str_len_max),*) "+line[2:vali-1].strip()+"\nexit\n")
				f90sbr2.write(word)
			elif line[0:2] == 'f_' and vali>0:
				word = ("!**real**\n\t"
				+"call getValue(ifound,iflag,\'"+section+"\',\'"+line[:vali-1].strip()+"\',"
				+"iniStr,\'"+line[vali:].strip()+"\')\n\t"
				+conditional_pre+"read(iniStr,*) "+line[2:vali-1].strip()+"\n\n\t")
				f90sbr.write(word)
				word = ("!**real**\n"
				+"case(\'"+line[:vali-1].strip()+"=\',\'"+line[2:vali-1].strip()+"=\')\n"
				+"read(name(istr+1:str_len_max),*) "+line[2:vali-1].strip()+"\nexit\n")
				f90sbr2.write(word)
			elif line[0:2] == 'i_' and vali>0:
				word = ("!**integer**\n\t"
				+"call getValue(ifound,iflag,\'"+section+"\',\'"+line[:vali-1].strip()+"\',"
				+"iniStr,\'"+line[vali:].strip()+"\')\n\t"
				+conditional_pre+"read(iniStr,*) "+line[2:vali-1].strip()+"\n\n\t")
				f90sbr.write(word)
				word = ("!**integer**\n"
				+"case(\'"+line[:vali-1].strip()+"=\',\'"+line[2:vali-1].strip()+"=\')\n"
				+"read(name(istr+1:str_len_max),*) "+line[2:vali-1].strip()+"\nexit\n")
				f90sbr2.write(word)
			elif line[0:2] == 's_' and vali>0:
				word = ("!**string**\n\t"
				+"call getValue(ifound,iflag,\'"+section+"\',\'"+line[:vali-1].strip()+"\',"
				+"iniStr,\'"+line[vali:].strip()+"\')\n\t"
				+conditional_pre+line[2:vali-1].strip()+"= trim(adjustl(iniStr))"+"\n\n\t")
				f90sbr.write(word)
				word = ("!**string**\n"
				+"case(\'"+line[:vali-1].strip()+"=\',\'"+line[2:vali-1].strip()+"=\')\n"
				+line[2:vali-1].strip()+" = name(istr+1:str_len_max)"+"\nexit\n")
				f90sbr2.write(word)
			elif line[0] != ';' and vali>0:
				word = ("call getValue(ifound,iflag,\'"+section+"\',\'"+line[:vali-1].strip()+"\',"
				+"iniStr,\'"+line[vali:].strip()+"\')\n\t"
				+conditional_pre+"read(iniStr,*) "+line[:vali-1].strip()+"\n\n\t")
				f90sbr.write(word)
				word = ("case(\'"+line[:vali-1].strip()+"=\',\'"+line[2:vali-1].strip()+"=\')\n"
				+"read(name(istr+1:str_len_max),*) "+line[2:vali-1].strip()+"\n\texit\n")
				f90sbr2.write(word)

f90sbr2.write("\ncase default\n\nend select\nenddo")
f90sbr2.write("\nend subroutine readcmdinput")
f90sbr.write("\nend subroutine readini")
f90sbr.close()
f90sbr2.close()