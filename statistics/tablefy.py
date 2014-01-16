#!/usr/bin/env python

from sys import argv

def parse(section):
	lines = section.split("\n")
	vars = {}
	for l in lines[1:]:
		values = l.split(":")
		if len(values) < 2:
			continue
		vars[values[0]] = float(values[1])
	return vars

def latexify(table):
	s = ""
	for t in table:
		s += t + " & "
		for c in table[t][:-1]:
			s += "%.3f"%(c) + " & "
		s += "%.3f"%(table[t][-1]) + " \\\\\n"
	return s

tables = [{},{}]
for a in argv[1:]:
	fl = open(a,"rt")
	text = fl.read()
	fl.close()
	
	sections = text.split("######")
	for i in range(len(sections)):
		s = sections[i]
		parsed = parse(s)
		for p in parsed:
			if p in tables[i]:
				tables[i][p].append(parsed[p])
			else:
				tables[i][p] = [parsed[p]]
	
# delta and percentual delta
for t in tables:
	for c in t:
		a = t[c][0]
		b = t[c][1]
		t[c].append(a-b)
		t[c].append(100*(a-b)/b)
		
for m in map(latexify,tables):
	print(m)