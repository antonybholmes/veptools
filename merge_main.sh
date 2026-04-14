msg=$1 #"Bug fixes and updates."
type="fix"
branch="dev"

if [[ -z "${msg}" ]]
then
	msg="Bug fixes and updates."
fi


OPTSTRING="t:m:b:"

while getopts ${OPTSTRING} opt
do
	case ${opt} in
  	t)
    	type=$OPTARG
      	;;
	m)
    	msg=$OPTARG
      	;;
	b)
      	branch=$OPTARG
      	;;
    ?)
      echo "Invalid option: -${OPTARG}."
      exit 1
      ;;
  esac
done

echo "${type}: ${msg}"
echo ${branch}

./commit.sh -t "${type}" -m "${msg}" -b dev

git switch main
git merge dev -m "${type}: ${msg}"

#git push -u origin main
./commit.sh -t "${type}" -m "${msg}" -b main

git switch dev
