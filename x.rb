lines = File.open("primates.html", "r").readlines
lines.each do |line|
  line.chomp!
  if (line.index("c:/genome_re") == nil) then
    print "#{line}\n"
  else
    data = line.split('"')
    stuff = data[1].split("/")
    print "#{data[0]}\"./images/#{stuff[6]}\"#{data[2]}\n"
  end
end
