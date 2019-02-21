# 总结perl产生随机数的rand和srand函数
## [引用](http://blog.sina.com.cn/s/blog_679686370102wowz.html)

srand(expr1)函数用于根据expr产生一个种子，此种子用于rand函数产生随机数。rand(expr2)函数用于根据种子产生一个随机数。由此可知，srand函数必定先于rand函数运行。

一、如果srand函数产生的种子是确定的，那么rand函数产生的随机数的序列就是确定的。像下面程序：

srand(0);

foreach (1..20)

{

        $rand_franction = rand(100);#产生一个介于[0,100)的随机分数。

        $rand_integer = int(rand(100));#产生一个介于[0,100)的整数。

}

无论上面的程序运行多少遍，每一遍生成的随机数序列都是相同的。

二、如果srand函数产生的种子是不确定的，那么rand函数产生的随机数序列也是不确定的。像下面程序：

srand();#或者利用系统时间产生种子：srand(time());

foreach (1..20)

{

        $rand_franction = rand(100);#产生一个介于[0,100)的随机分数。

        $rand_integer = int(rand(100));#产生一个介于[0,100)的整数。

}

上面的程序每次运行产生的随机数的序列都是不同的；当然如果碰巧某两次运行产生的种子是相同的，那么随机数序列也是相同的，这种可能性是存在的。

三、不显式地调用srand函数产生种子时，在第一次运行rand函数前，srand函数会被自动调用，生成一个随机的种子后，才执行rand函数产生随机数。这种情况和显式地执行srand的效果是一样的。如下面程序：

#srand();

foreach (1..20)

{

        $rand_franction = rand(100);#产生一个介于[0,100)的随机分数。

        $rand_integer = int(rand(100));#产生一个介于[0,100)的整数。

}

四、无论是显式地调用srand，还是利用rand去自动调用srand函数，srand函数可以只运行一次，而不管rand函数被执行多少次。
