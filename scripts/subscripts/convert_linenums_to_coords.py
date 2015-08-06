import sys

uniqs = sys.stdin
linenums = open(sys.argv[1])

pos = -1

for line in linenums:

	(i, l) = line.strip().split('\t')
	(i, l) = (int(i), int(l))

	# if the line number is greater than the current
	# uniquely mapping segment, get the next one.

	while i > pos:

		rline = uniqs.readline()
		if not rline:
			break

		(contig, start, end) = rline.strip().split('\t')
		(start, end) = (int(start), int(end))

		pos += end - start

	# sample a second tag from the fragment length distribution

	s = end - (pos - i + 1)

	print "%s\t%d\t%d" % (contig, s, s + 1)
	print "%s\t%d\t%d" % (contig, s + l, s + l + 1);

# If we don't read through the rest of what is being piped in, that part of the pipe can end up
# with a 141 status in bash and cause things to exit due to error
for line in uniqs:
    pass

linenums.close()
sys.exit(0)
