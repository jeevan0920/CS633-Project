echo "Wrtie active ip's to the hostfile Loading ... "
rm hostfile 1>/dev/null 2>/dev/null
ifconfig | grep inet | head -1 | awk -F' ' '{print $2}' > hostfile
n=4
m=$(ifconfig | grep inet | head -1 | awk -F' ' '{print $2}' | awk -F '.' '{print $4}')
while [ $n -le 30 ]; do
        node_ip=172.27.19.$n

        #do not add  current  ip
        if [ $n -eq $m ];then
                n=`expr $n + 1`
                continue
        fi

        #add the IP to host file
        ping -c 1 $node_ip -W 1 > /dev/null
        if [ $? -eq 0 ];then
                echo $node_ip:4 >> hostfile
        fi
        n=`expr $n + 1`
done

echo "Hostfiles are created "
