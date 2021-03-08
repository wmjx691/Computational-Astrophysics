! 驚嘆號之後是註解
program main ! 這行可以省略，但是寫大程式的時候會發生混亂
    write (*,*) "hello, world!" ! 第一個* 表示輸出縮排使用內定值，第二個* 表示不指定輸出格式
    write (unit = *, fmt = * ) "hello, world!" ! 做和上一行一樣的事
    stop ! 這行代表程式結束，可以省略
end program main ! end之後的program main也可以省略，但寫上是比較嚴謹