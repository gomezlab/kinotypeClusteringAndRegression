if command -v npm > /dev/null
then
    sleep 1
else
    echo "Error. Please install npm (likely with sudo apt get install npm)"
    exit 1
fi

# be sure we've installed http-server
sudo npm install -g http-server

# run the packing script with the given input file
python pack.py $1

# check if runtime error
ret=$?
if [ $ret -ne 0 ]; then
    echo "Parsing error. Is the file correctly formatted as tab-separated?"
    exit 1
fi

# run on port 8080
npm start